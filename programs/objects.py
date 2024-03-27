from datetime import datetime
from math import sqrt, tan
import numpy as np

from config import *
from mechanics import Mechanics, Earth
from MathPy import Math, Triangulate
from tools import Parser, Filer, Grapher


class Satellite:
    '''
    Искусственный спутник

    Параметры:
    ------
        `type`: `str` - тип спутника
    '''
    def __init__(self, type: str, t0=0, parameters=parameters):
        self.type = type
        # Elements
        self.a = parameters[type]['a']
        self.e = parameters[type]['e']
        self.i = parameters[type]['i']
        self.w = parameters[type]['w']
        self.Omega = parameters[type]['Omega']
        self.M0 = parameters[type]['M0']
        self.T = parameters[type]['T']

        self.t0 = t0
        self.n = sqrt(Mechanics.Mu / self.a ** 3)
        # Initial coords
        self.orbital_coords, self.orbital_velocities = Mechanics.get_orbital_coords(self, t0)
        self.coords, self.velocities = Mechanics.get_coords(self, t0)

    def __repr__(self):
        return repr(f'<Satellite: {self.type}>')

    def evolution(self):
        '''Эволюция динамики за период обращения'''
        dt = 0.01 * self.T # c
        t = self.t0
        while t <= self.t0 + self.T:
            self.coords, self.velocities = Mechanics.get_coords(self, t)
            x, y, z = self.coords
            yield np.array([x, y, z])
            t += dt

    def route(self, date: str):
        '''
        Трасса спутника
        '''
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
            coords.append(Math.get_lmd_phi(y))
            t += dt
        return coords
    

class Point:
    def __init__(self, coords: tuple):
        self.x, self.y, self.z = coords
        self.coords = coords


class Object:
    '''
    Объект, воздействующий на спутник
    '''
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
    fm = 332946.04877304828 * Earth.fm
    

class Moon(Object):
    fm = Earth.fm / 81.30056822149722
    

class GLONASS:
    '''
    Система глобального позиционирования
    '''
    def __init__(self):
        self._parser = Parser()
        self.data = self._parser.parse_data()

    def current_position(self, **kwargs):
        '''
        Положение спутников в настоящий момент времени (по данным https://glonass-iac.ru/glonass/)
        '''
        grapher = Grapher(custom_rcParams, None)
        
        result = self._parser.get_current_position()
        lmd, phi = zip(*result.values())
        text = list(result.keys())

        grapher.ax.scatter(lmd, phi, **marker)
        grapher.ax.set_xlabel('$\lambda, °$')
        grapher.ax.set_ylabel('$\phi, °$')
        
        for i in range(len(lmd)):
            grapher.ax.annotate(text[i], (lmd[i], phi[i] + 2))

        grapher.ax.set_title(f'{datetime.now()}')
        grapher.show()

    def __get_elems(self, H_omega: float, data: dict) -> tuple[float]:
        ecc = data['e']
        Omega = data['Lomega'].rad + H_omega

        v_omega = -data['W']
        E_omega = 2 * np.arctan2(tan(v_omega/2) * sqrt((1 - ecc)/(1 + ecc)), 1)
        M_omega = E_omega - ecc * np.sin(E_omega)
        a = (data['Tapp'] / (2 * np.pi) * sqrt(Mechanics.Mu)) ** (0.66)

        return(a, Omega, M_omega)

    def get_sat_elements(self):
        '''
        Определение элементов орбит спутников 
        '''
        parameters = {}
        for satellite, data in self.data.items():
            date = f'{data["datetime"]} 00:00:00'
            JD = Math.get_JD(date) - 3/24
            JD_omega = JD + data['Tomega']/86400
            H_omega = Math.sid2000(JD_omega)

            a, Omega, M_omega = self.__get_elems(H_omega, data)
            parameters[satellite] = {
                'a': a,
                'e': data['e'],
                'i': data['i'],
                'T': data['Tapp'],
                'Omega': Omega,
                'w': data['W'],
                'M_omega': M_omega,
                'JD_omega': JD_omega,
                'date': date,
            }
        return parameters

    def get_satellites(self, time: str) -> list[Satellite]:
        JD0 = Math.get_JD(time)
        parameters = self.get_sat_elements()

        sats = []
        for sat in parameters:
            JD_omega = parameters[sat]['JD_omega']

            dt = (JD0 - JD_omega) * 86400
            n = 2 * np.pi / parameters[sat]['T']
            M = parameters[sat]['M_omega'] + n * dt
            parameters[sat]['M0'] = M
            satellite = Satellite(type=sat, parameters=parameters)
            sats.append(satellite)
        return sats
            
    def get_positions(self, time: str) -> list[Satellite]:
        sats = self.get_satellites(time)
        grapher = Grapher(custom_rcParams, projection=None)

        for sat in sats:
            geo_coords = Math.get_lmd_phi(sat.coords)
            grapher.print(geo_coords, c='red')
            grapher.ax.set_title(time)
            grapher.ax.annotate(int(sat.type), (geo_coords[0].decimal, 
                                           geo_coords[1].decimal + 2))
        grapher.show() 

    def get_routes(self, time: str) -> list[Satellite]:
        sats = self.get_satellites(time)
        grapher = Grapher(custom_rcParams, projection=None)

        for sat in sats:
            route = sat.route(time)
            initial_position = route[0]
            
            grapher.print(initial_position, c='red')
            grapher.ax.annotate(int(sat.type), (initial_position[0].decimal, 
                                           initial_position[1].decimal + 2))
            grapher.print(route, s=1, c='k')
        grapher.show() 
    
    def get_orbits(self, time: str):
        sats = self.get_satellites(time)
        earth = Earth(EARTH_PATH)
        grapher = Grapher(custom_rcParams, projection='3d')
        for sat in sats:
            evolution = np.array(list(sat.evolution())).T
            coords = np.array(list(sat.coords)).T
            
            grapher.print(evolution, c='black')
            grapher.print(earth.coords.T, c='royalblue')
            grapher.print(coords, c='red', s=10)

        grapher.show()


# g = GLONASS()
# g.current_position()
# g.get_positions('14.03.2024 11:58:04')
# g.get_routes('14.03.2024 11:58:04')
# g.get_orbits('14.03.2024 10:16:34')