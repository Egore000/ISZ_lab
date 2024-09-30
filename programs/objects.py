from datetime import datetime
from math import sqrt, tan, sin, cos
import numpy as np
from loguru import logger

from config import *
from mechanics import Mechanics, Earth
from MathPy import Math, Triangulate
from tools import Parser, Filer, Grapher


class Satellite:
    """
    Искусственный спутник

    Параметры:
    ------
        `type`: `str` - тип спутника

        `t0`: `float` - начальная эпоха (юлианская дата)
        
        `parameters`: `dict` - параметры спутника
    """

    def __init__(self, type: str, t0: float = 0, parameters: dict = parameters):
        self.type = type
        # Elements
        self.a = parameters[type]["a"]
        self.e = parameters[type]["e"]
        self.i = parameters[type]["i"]
        self.w = parameters[type]["w"]
        self.Omega = parameters[type]["Omega"]
        self.M0 = parameters[type]["M0"]
        self.T = parameters[type]["T"]

        self.t0 = t0
        self.n = sqrt(Mechanics.Mu / self.a ** 3)
        # Initial coords
        self.orbital_coords, self.orbital_velocities = Mechanics.get_orbital_coords(
            self, t0)
        self.coords, self.velocities = Mechanics.get_coords(self, t0)

        self.x, self.y, self.z = self.coords
        self.Vx, self.Vy, self.Vz = self.velocities

    def __repr__(self):
        return repr(f"<Satellite: {self.type}>")

    def evolution(self):
        """Эволюция динамики за период обращения"""
        dt = 0.01 * self.T  # c
        t = self.t0
        while t <= self.t0 + self.T:
            self.coords, self.velocities = Mechanics.get_coords(self, t)
            x, y, z = self.coords
            yield np.array([x, y, z])
            t += dt

    def route(self, date: str) -> list[tuple[Angles, Angles]]:
        """Трасса спутника"""
        
        jd = Math.get_JD(date)      # Юлианская дата
        (_, _, t0) = Math.get_daytime(date)
        t0 -= int(t0)
        t0 *= 86400

        dt = 0.01 * self.T
        t = 0
        coords = []
        while t < self.T:
            H0 = Math.sid2000(jd + t/86400)
            H = H0 + Mechanics.w * t0

            y = Mechanics.TRS(self, H, t)
            # coords.append(Math.get_lmd_phi(y))

            lmd, phi = Math.get_lmd_phi(y)
            lmd -= Angles(180)

            if lmd < Angles(-180):
                lmd += Angles(360)

            coords.append((lmd, phi))

            t += dt
        return coords


class Point:
    """Точка в пространстве"""
    def __init__(self, coords: tuple):
        self.x, self.y, self.z = coords
        self.coords = coords


class Object:
    """Объект, воздействующий на спутник"""

    def __init__(self, coords: tuple):
        self.coords = coords

    def disturbing_function(self, satellite: Satellite):
        (x, y, z) = satellite.coords
        (x_, y_, z_) = self.coords

        dx = (x-x_, y-y_, z-z_)
        r_ = Math.radius(self.coords)
        delta = Math.radius(dx)

        return self.fm * (1 / delta - (x*x_ + y*y_ + z*z_) / r_**3)


class Sun(Object):
    """Солнце"""
    fm = 332946.04877304828 * Earth.fm


class Moon(Object):
    """Луна"""
    fm = Earth.fm / 81.30056822149722


class GLONASS:
    """
    Система глобального позиционирования ГЛОНАСС
    """

    def __init__(self):
        self._parser = Parser()
        self.data = self._parser.parse_data

    def current_position(self, **kwargs):
        """Положение спутников в настоящий момент времени
          (по данным https://glonass-iac.ru/glonass/)"""
        
        grapher = Grapher(custom_rcParams, None)

        result = self._parser.get_current_position()
        lmd, phi = zip(*result.values())
        text = list(result.keys())

        grapher.ax.scatter(lmd, phi, **marker)
        grapher.ax.set_xlabel("$\lambda, °$")
        grapher.ax.set_ylabel("$\phi, °$")

        for i in range(len(lmd)):
            grapher.ax.annotate(text[i], (lmd[i], phi[i] + 2))

        grapher.ax.set_title(f"{datetime.now()}")
        grapher.show()

    def __get_elements(self, H_omega: float, data: dict) -> tuple[float]:
        ecc = data["e"]
        Omega = data["Lomega"].rad + H_omega

        v_omega = -data["W"]
        E_omega = 2 * np.arctan2(tan(v_omega/2) * sqrt((1 - ecc)/(1 + ecc)), 1)
        M_omega = E_omega - ecc * np.sin(E_omega)
        a = (data["Tapp"] / (2 * np.pi) * sqrt(Mechanics.Mu)) ** (0.66)

        return (a, Omega, M_omega)

    def get_sat_elements(self):
        """Определение элементов орбит спутников"""
        
        parameters = {}
        for satellite, data in self.data().items():
            date = f"{data["datetime"]} 00:00:00"
            JD = Math.get_JD(date) - 3/24
            JD_omega = JD + data["Tomega"]/86400
            H_omega = Math.sid2000(JD_omega)

            a, Omega, M_omega = self.__get_elements(H_omega, data)
            parameters[satellite] = {
                "a": a,
                "e": data["e"],
                "i": data["i"],
                "T": data["Tapp"],
                "Omega": Omega,
                "w": data["W"],
                "M_omega": M_omega,
                "JD_omega": JD_omega,
                "date": date,
            }
        return parameters

    def get_satellites(self, time: str):
        """Получение спуников группировки ГЛОНАСС"""

        JD0 = Math.get_JD(time)
        parameters = self.get_sat_elements()

        for sat in parameters:
            JD_omega = parameters[sat]["JD_omega"]

            dt = (JD0 - JD_omega) * 86400
            n = 2 * np.pi / parameters[sat]["T"]
            M = parameters[sat]["M_omega"] + n * dt
            parameters[sat]["M0"] = M
            satellite = Satellite(type=sat, parameters=parameters)
            yield satellite

    def get_positions(self, time: str):
        """Построение графика с положением спутников ГЛОНАСС в заданный момент времени"""

        satellites = self.get_satellites(time)
        grapher = Grapher(custom_rcParams, projection=None)

        for sat in satellites:
            geo_coords = Math.get_lmd_phi(sat.coords)
            grapher.print(geo_coords, c="red")
            grapher.ax.annotate(int(sat.type), (geo_coords[0].decimal,
                                                geo_coords[1].decimal + 2))
        grapher.ax.set_title(time)
        grapher.show()

    def get_routes(self, time: str) -> list[Satellite]:
        """Построение трасс спутников ГЛОНАСС"""

        satellites = self.get_satellites(time)
        grapher = Grapher(custom_rcParams, projection=None)

        for sat in satellites:
            route = sat.route(time)
            initial_position = route[0]

            grapher.ax.annotate(int(sat.type), (initial_position[0].decimal,
                                                initial_position[1].decimal + 2))
            grapher.print(route, s=1, c="k")
            grapher.print(initial_position, c="red")

        grapher.ax.set_title(time)
        grapher.show()

    def get_orbits(self, time: str):
        """Построение орбит спутников ГЛОНАСС в 3D проекции"""

        satellite = Satellite(type="Геостационарный")

        jd = Math.get_JD(time)
        t = jd * 86400

        satellite.coords, _ = Mechanics.get_coords(satellite, t)

        satellites = self.get_satellites(time)
        earth = Earth(EARTH_PATH)
        grapher = Grapher(custom_rcParams, projection="3d")

        grapher.ax.scatter(satellite.coords[0], satellite.coords[1], c="red")
        for sat in satellites:
            evolution = np.array(list(sat.evolution())).T
            coords = np.array(list(sat.coords)).T

            grapher.print(evolution, c="black")
            grapher.print(coords, c="red", s=10)

        grapher.print(earth.coords.T, c="royalblue")
        grapher.ax.set_title(time)
        grapher.show()

    def _set_params(self, grapher, title):
        grapher.ax.set_title(title)
        grapher.ax.set_xlim(-50000, 50000)
        grapher.ax.set_ylim(-50000, 50000)
        grapher.ax.set_xlabel("$\it{y}$, км")
        grapher.ax.set_ylabel("$\it{z}$, км")
        grapher.ax.set_aspect("equal", adjustable="box")

    def visibility_areas(self, time: str):
        """Построение зон видимости спутников ГЛОНАСС"""
        satellite = Satellite(type="Геостационарный")

        jd = Math.get_JD(time)
        t = jd * 86400

        satellite.coords, _ = Mechanics.get_coords(satellite, t)

        satellites = self.get_satellites(time)
        earth = Earth(EARTH_PATH)
        grapher = Grapher(custom_rcParams, projection=None)

        grapher.print(earth.coords[:, 1:], c="royalblue")
        for sat in satellites:
            x, y, z = sat.coords

            grapher.print((y, z), c="black")
            grapher.ax.annotate(int(sat.type), (y, z + 2))

        grapher.print(satellite.coords[1:], c="red")

        self._set_params(grapher, time)
        grapher.show()


class Observer:
    """Наблюдатель"""
    def __init__(self, lat: Angles, lon: Angles, height: float):
        self._lat = lat
        self._lon = lon
        self._h = height

    @property
    def coords(self) -> tuple[Angles, Angles]:
        """Сферические координаты наблюдателя"""
        return (self._lon, self._lat)
    
    def coords_CRS(self, sidereal_time: Angles) -> tuple[float, float, float]:
        """Координаты наблюдателя в небесной системе отсчёта"""
        return self._get_coords_CRS(sidereal_time)

    def local_sidereal_time(self, sidereal_time: Angles) -> Angles:
        """Локальное звёздное время"""
        return self._get_local_sidereal_time(sidereal_time)

    def _get_coords_CRS(self, sidereal_time: Angles) -> tuple[float, float, float]:
        H = self._get_local_sidereal_time(sidereal_time)

        x = (Earth.Radius + self._h) * cos(self._lat) * cos(H)
        y = (Earth.Radius + self._h) * cos(self._lat) * sin(H)
        z = (Earth.Radius + self._h) * sin(self._lat)
        
        return (x, y, z)
        
    def _get_local_sidereal_time(self, sidereal_time: Angles) -> Angles:
        return sidereal_time + self._lon