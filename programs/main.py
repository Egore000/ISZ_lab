from math import sqrt

from config import *
from tools import *
from mechanics import *


class Satellite:
    def __init__(self, t=T0):
        # Elements
        self.a = a
        self.e = e
        self.i = i
        self.w = w
        self.Omega = Omega
        self.M0 = M0
        self.t0 = T0
        self.T = T * 3600
        self.n = sqrt(Mechanics.Mu / self.a ** 3)
        # Coords
        self.orbital_coords, self.orbital_velocities = Mechanics.get_orbital_coords(self, t)
        self.coords, self.velocities = Mechanics.get_coords(self, t)

    def evolution(self):
        '''Эволюция динамики за период обращения'''
        t = self.t0
        dt = 0.01 * self.T # c
        while t <= self.t0 + self.T:
            yield (t, Mechanics.get_coords(self, t)[0])
            t += dt

    # def route(self, date: str):
    #     jd = Math.get_JD(date)      # Юлианская дата
    #     (_, _, t0) = Math.get_daytime(date)    
    #     t0 -= int(t0) 
    #     H0 = Math.sid2000(jd)       # Звёздное время для данной юлианской даты
    #     t = self.t0
    #     dt = 0.01 * self.T
    #     t0 *= 86400

    #     while t < self.t0 + self.T:
    #         t0 += dt           
    #         print(t0)
    #         H = H0.rad + Mechanics.w * t0

    #         x, v = Mechanics.get_coords(self, t)
    #         x = np.array(x)
    #         y = Mechanics.transition(H, x)
    #         yield Math.get_lmd_phi(y)
    #         t += dt
            
    def route(self, date: str):
        jd = Math.get_JD(date)                              # Юлианская дата 
        jd00 = Math.get_JD(date.split()[0] + ' 00:00')      # Юлианская дата на 00:00
        jd2000 = Math.get_JD('01.01.2000 00:00')            # Эпоха J2000.0
        
        t0 = (jd - jd00) * 86400

        dt = 0.01 * self.T
        t = 0
        while t < self.T:               
            d = (jd + t / 86400) - jd2000
            Tu = d / 36525

            H0 = 24110.54841 + 8640184.812866 * Tu + 0.093104 * Tu**2 - 6.21e-6 * Tu**3
            H0 = H0 / 3600 * pi/180

            H = H0 + Mechanics.w * t0
            H = Angles(rad=H).reduce()
  
            x, v = Mechanics.get_coords(self, t)
            x = np.array(x)
            y = Mechanics.transition(H, x)
            yield Math.get_lmd_phi(y)
            t += dt


def orbit():
    grapher = Grapher(custom_rcParams)
    earth = Earth(EARTH_PATH)
    satellite = Satellite()

    time, evolution = zip(*satellite.evolution())
    grapher.print(evolution, c='black')
    grapher.print(earth.coords)
    grapher.print(satellite.coords, c='red', s=10)

    grapher.show()
    

def route():
    satellite = Satellite()

    route = satellite.route(date)
    lmd, phi = zip(*route)
    lmd = list(map(lambda x: x.decimal, lmd))
    phi = list(map(lambda x: x.decimal, phi))

    plt.scatter(lmd, phi, s=1)
    plt.show()


if __name__=="__main__":
    route()
