from math import sin, cos, sqrt, pi
import numpy as np
import numpy.typing as npt

from AnglesPy import Angles
from MathPy import Math
from config import *
from tools import Filer


class Earth:
    Mass = 2e24
    Radius = 6378.1363
    fm = 398600.5

    def __init__(self, file: str):
        self.coords = Filer.read(file)


class Mechanics:
    """
    Класс для работы с формулами из небесной механики
    """

    Mu = 3.986004418E5
    f = 6.67e-17            # H * km**2 / kg**2
    w = 7.2922115e-5        # скорость среднего звёздного вращения Земли
    ae = 149_597_870.7      # km
    c = 299_792.458         # km/s

    @staticmethod
    def rotation_matrix(axis: str, angle: Angles | float) -> npt.NDArray:
        """Матрица поворота"""
        if axis in ("x", "Ox"):
            return np.array([
                [1, 0, 0],
                [0, cos(angle), -sin(angle)],
                [0, sin(angle), cos(angle)]
            ])

        elif axis in ("y", "Oy"):
            return np.array([
                [cos(angle), 0, sin(angle)],
                [0, 1, 0],
                [-sin(angle), 0, cos(angle)]
            ])

        elif axis in ("z", "Oz"):
            return np.array([
                [cos(angle), sin(angle), 0],
                [-sin(angle), cos(angle), 0],
                [0, 0, 1]
            ])

    def get_elements(self, coords: list, velocities: list):
        """
        Переход от координат и скоростей к элементам орбиты

        Параметры::
        ------
            `coords` - массив координат

            `velocities` - массив скоростей

        Выходные значения::
        ------
        `(ecc, i, a, Omega, w, M)`
        """
        x, y, z = coords
        Vx, Vy, Vz = velocities

        r = Math.radius(coords)
        V2 = Math.radius(velocities)
        h = V2 / 2 - self.Mu / r

        c1 = y * Vz - z * Vy
        c2 = z * Vx - x * Vz
        c3 = x * Vy - y * Vx

        l1 = -self.Mu * x/r + Vy * c3 - Vz * c2
        l2 = -self.Mu * y/r + Vz * c1 - Vx * c3
        l3 = -self.Mu * z/r + Vx * c2 - Vy * c1

        c = Math.radius((c1, c2, c3))
        l = Math.radius((l1, l2, l3))

        a = -self.Mu / (2 * h)
        ecc = l / self.Mu
        i = np.arctan2(sqrt(1 - c3**2 / c**2), c3 / c)

        if i == 0:
            i += 1e-12

        Omega = np.arctan2(c1/(c*sin(i)), -c2/(c*sin(i)))
        w = np.arctan2(l3/(l*sin(i)), l1/l*cos(Omega) + l2/l*sin(Omega))
        E0 = np.arctan2((x*Vx + y*Vy + z*Vz) /
                        (ecc*sqrt(self.Mu*a)), (1-r/a)/ecc)
        M = E0 - ecc*sin(E0)

        i, Omega, w, M = list(map(lambda x: Angles(rad=x), [i, Omega, w, M]))
        return (ecc, i, a, Omega, w, M)

    @staticmethod
    def __get_parameters(satellite):
        """Получение параметров"""

        a1 = cos(satellite.w)*cos(satellite.Omega) - \
            sin(satellite.w)*sin(satellite.Omega)*cos(satellite.i)
        b1 = cos(satellite.w)*sin(satellite.Omega) + \
            sin(satellite.w)*cos(satellite.Omega)*cos(satellite.i)
        c1 = sin(satellite.w)*sin(satellite.i)

        a2 = -sin(satellite.w)*cos(satellite.Omega) - \
            cos(satellite.w)*sin(satellite.Omega)*cos(satellite.i)
        b2 = -sin(satellite.w)*sin(satellite.Omega) + \
            cos(satellite.w)*cos(satellite.Omega)*cos(satellite.i)
        c2 = cos(satellite.w)*sin(satellite.i)

        return (a1, b1, c1), (a2, b2, c2)

    @staticmethod
    def anomaly(satellite, t: float):
        """Вычисление эксцентрической аномалии спутника `satellite`"""

        M = satellite.n * (t - satellite.t0) + satellite.M0
        dif = 1
        E0 = M
        while abs(dif) > 1e-13:
            E = M + satellite.e * sin(E0)
            dif = E - E0
            E0 = E
        return E

    @staticmethod
    def get_orbital_coords(satellite, t: float):
        """Вычисление орбитальных координат и скоростей спутника `satellite`"""

        E = Mechanics.anomaly(satellite, t)

        xi = satellite.a * (cos(E) - satellite.e)
        eta = satellite.a * sqrt(1 - satellite.e**2) * sin(E)

        Vxi = -satellite.a * satellite.n * sin(E)/(1 - satellite.e * cos(E))
        Veta = satellite.a * satellite.n * \
            sqrt(1 - satellite.e**2) * cos(E)/(1 - satellite.e * cos(E))

        return (xi, eta), (Vxi, Veta)

    @staticmethod
    def get_coords(satellite, t: float):
        """Вычисление координат спутника `satellite`"""

        (a1, b1, c1), (a2, b2, c2) = Mechanics.__get_parameters(satellite)

        (xi, eta), (Vxi, Veta) = Mechanics.get_orbital_coords(satellite, t)

        x = a1*xi + a2*eta
        y = b1*xi + b2*eta
        z = c1*xi + c2*eta

        Vx = a1*Vxi + a2*Veta
        Vy = b1*Vxi + b2*Veta
        Vz = c1*Vxi + c2*Veta

        return (x, y, z), (Vx, Vy, Vz)

    @staticmethod
    def transition(h: float, x):
        """Переход во вращающуюся систему координат"""
        A = Mechanics.rotation_matrix("z", h)
        return A @ x

    @staticmethod
    def geopotential(coords: tuple[float],
                     n1: int, n2: int) -> float:
        """
        Вычисление геопотенциала в заданной точке `coords` с учётом
        влияния гармоник геопотенциала 

        `n1` - начальная гармоника

        `n2` - конечная гармоника
        """
        r0 = Earth.Radius
        lmd, phi = Math.get_lmd_phi(coords)
        r = Math.radius(coords)

        coefs = Filer.read_coefficients(COEFS)

        coef = Earth.fm / r

        sum = 0
        for n in range(n1, n2 + 1):
            for m in range(n + 1):
                (Cnm, Snm) = coefs[(n, m)]
                sum += (r0/r) ** (n + 1) * Math.Norm_LPnm(n, m, sin(phi))\
                    * (Cnm * cos(m * lmd) + Snm * sin(m * lmd))
        return coef * sum

    @staticmethod
    def TRS(satellite, H: float, t: float):
        """Координаты спутника во вращающейся системе координат"""
        x, v = Mechanics.get_coords(satellite, t)
        x = np.array(x)
        y = Mechanics.transition(H, x)
        return y

    @staticmethod
    def HTS(vector: npt.NDArray, local_sidereal_time: Angles, lat: Angles) -> npt.NDArray:
        H = local_sidereal_time

        A = Mechanics.rotation_matrix("z", H)
        M = Mechanics.rotation_matrix("y", Angles(90) - lat)

        return M @ A @ vector

    @staticmethod
    def from_HTS_to_CRS(vector: npt.NDArray, local_sidereal_time: Angles, lat: Angles) -> npt.NDArray:
        H = local_sidereal_time

        A = Mechanics.rotation_matrix("z", H)
        M = Mechanics.rotation_matrix("y", Angles(90) - lat)

        return np.linalg.inv(A) @ np.linalg.inv(M) @ vector
    
    @staticmethod
    def from_CRS_to_TRS(vector: npt.NDArray, time: str) -> tuple[Angles, Angles]:
        jd = Math.get_JD(time) 
        (_, _, t0) = Math.get_daytime(time)
        t0 -= int(t0)
        t0 *= 86400     # Юлианская дата

        H0 = Math.sid2000(jd)
        H = Angles(rad=H0 + Mechanics.w * t0)

        trs = Mechanics.transition(H, vector)
        return Math.get_lmd_phi(trs)


