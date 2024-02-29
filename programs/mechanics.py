from math import sin, cos, sqrt, pi
import numpy as np

from AnglesPy import Angles
from MathPy import Math
from config import * 
from tools import Filer, Earth

class Mechanics:
    '''Класс для работы с формулами из небесной механики'''
    Mu = 3.986004418E5
    f = 6.67e-17            # H * km**2 / kg**2
    w = 7.2922115e-5        # скорость среднего звёздного вращения Земли
    ae = 149_597_870.7        # km
    def get_elements(self, coords: list, velocities: list):
        '''
        Переход от координат и скоростей к элементам орбиты

        Параметры::
        ------
            `coords` - массив координат
        
            `velocities` - массив скоростей

        Выходные значения::
        ------
        `(ecc, i, a, Omega, w, M)`
        '''
        x, y, z = coords[0], coords[1], coords[2]
        Vx, Vy, Vz = velocities[0], velocities[1], velocities[2]

        r = sqrt(x**2 + y**2 + z**2)
        V2 = Vx**2 + Vy**2 + Vz**2
        h = V2/2 - self.Mu/r
        
        c1 = y*Vz - z*Vy
        c2 = z*Vx - x*Vz 
        c3 = x*Vy - y*Vx 

        l1 = -self.Mu*x/r + Vy*c3 - Vz*c2
        l2 = -self.Mu*y/r + Vz*c1 - Vx*c3
        l3 = -self.Mu*z/r + Vx*c2 - Vy*c1

        c = sqrt(c1**2 + c2**2 + c3**2)
        l = sqrt(l1**2 + l2**2 + l3**2)

        a = -self.Mu/(2*h)
        ecc = l/self.Mu
        i = np.arctan2(sqrt(1 - c3**2/c**2), c3/c)

        if i == 0:
            i += 1e-12

        Omega = np.arctan2(c1/(c*sin(i)), -c2/(c*sin(i)))
        w = np.arctan2(l3/(l*sin(i)), l1/l*cos(Omega) + l2/l*sin(Omega))
        E0 = np.arctan2((x*Vx + y*Vy + z*Vz)/(ecc*sqrt(self.Mu*a)), (1-r/a)/ecc)
        M = E0 - ecc*sin(E0)
        
        i, Omega, w, M = list(map(lambda x: Angles(rad=x), [i, Omega, w, M]))
        return (ecc, i, a, Omega, w, M)
    
    @staticmethod
    def __get_parameters(satellite):
        '''Получение параметров'''
        a1 = cos(satellite.w)*cos(satellite.Omega) - sin(satellite.w)*sin(satellite.Omega)*cos(satellite.i)
        b1 = cos(satellite.w)*sin(satellite.Omega) + sin(satellite.w)*cos(satellite.Omega)*cos(satellite.i)
        c1 = sin(satellite.w)*sin(satellite.i)

        a2 = -sin(satellite.w)*cos(satellite.Omega) - cos(satellite.w)*sin(satellite.Omega)*cos(satellite.i)
        b2 = -sin(satellite.w)*sin(satellite.Omega) + cos(satellite.w)*cos(satellite.Omega)*cos(satellite.i)
        c2 = cos(satellite.w)*sin(satellite.i)

        return (a1, b1, c1), (a2, b2, c2)

    @staticmethod
    def anomaly(satellite, t: float):
        '''Вычисление эксцентрической аномалии спутника `satellite`'''
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
        '''Вычисление орбитальных координат и скоростей спутника `satellite`'''
        E = Mechanics.anomaly(satellite, t)

        xi = satellite.a * (cos(E) - satellite.e)
        eta = satellite.a * sqrt(1 - satellite.e**2) * sin(E)

        Vxi = -satellite.a * satellite.n * sin(E)/(1 - satellite.e * cos(E))
        Veta = satellite.a * satellite.n * sqrt(1 - satellite.e**2) * cos(E)/(1 - satellite.e * cos(E))

        return (xi, eta), (Vxi, Veta)
    
    @staticmethod
    def get_coords(satellite, t: float):
        '''Вычисление координат спутника `satellite`'''
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
        '''Переход во вращающуюся систему координат'''
        A = np.array([
                    [cos(h), sin(h), 0], 
                    [-sin(h), cos(h), 0],
                    [0, 0, 1]
                    ])
        return A @ x
    
    @staticmethod
    def geopotential(coords: tuple[float], 
                     n1: int, n2: int) -> float:
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