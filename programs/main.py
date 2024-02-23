from math import sqrt

from config import *
from tools import *
from mechanics import *


class Satellite:
    '''
    Искусственный спутник

    Параметры:
    ------
        `type`: `str` - тип спутника
    '''
    def __init__(self, type: str, t0=0):
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
        coords = []
        while t <= self.t0 + self.T:
            self.coords, self.velocities = Mechanics.get_coords(self, t)
            coords.append((t, self.coords))
            t += dt
        return coords

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



def orbit(type: str):
    earth = Earth(EARTH_PATH)
    satellite = Satellite(type=type, t0=1e14)
    evolution = satellite.evolution()

    time, evolution = zip(*evolution)
    
    grapher = Grapher(custom_rcParams, projection='3d')
    grapher.print(evolution, c='black')
    grapher.print(earth.coords)
    print(satellite.coords)
    grapher.print(satellite.coords, c='red', s=10)
    grapher.fig.suptitle(type + ' спутник')
    grapher.show()
    

def route(type: str):
    satellite = Satellite(type=type)
    route = satellite.route(date)
    
    grapher = Grapher(custom_rcParams, projection=None)
    grapher.print(route, title=type + ' спутник', **marker)
    grapher.show()


if __name__=="__main__":
    # route('Тестовый')
    orbit('Тестовый')
    # route('Геостационарный')
    # orbit('Геостационарный')
