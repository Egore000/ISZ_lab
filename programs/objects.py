from math import sqrt
import numpy as np

from config import *
from mechanics import Mechanics, Earth
from MathPy import Math


class Satellite:
    '''
    Искусственный спутник

    Параметры:
    ------
        `type`: `str` - тип спутника
    '''
    def __init__(self, type: str, t0=0):
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

            x, v = Mechanics.get_coords(self, t)
            x = np.array(x)
            y = Mechanics.transition(H, x)
            coords.append(Math.get_lmd_phi(y))
            t += dt
        return coords


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
    

