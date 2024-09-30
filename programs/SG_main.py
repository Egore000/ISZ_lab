import numpy as np

from config import *
from mechanics import Mechanics
from MathPy import Math, Triangulate, Leveling
from objects import Satellite, Earth, Sun, Moon, Point


__doc__ = """
Задачи по спутниковой геодезии
"""

def disturb():
    """Определение возмущающей функции"""
    sun = Sun((
        -1.848432950986775E-01 * Mechanics.ae,
              8.860598024929441E-01 * Mechanics.ae,
              3.841027854355354E-01 * Mechanics.ae
              ))

    moon = Moon((
                (-1.863831425711169E-01 + 1.848432950986775E-01) * Mechanics.ae,
                (8.878791230277183E-01 - 8.860598024929441E-01) * Mechanics.ae,
                (3.850801888852833E-01 - 3.841027854355354E-01) * Mechanics.ae
                ))

    satellite = Satellite(type='Тестовый')
    satellite.coords = (-20710.0198816410957026,
                        13917.1445923690977837,
                        5402.62160300007461444)

    U = Earth.fm / Math.radius(satellite.coords)

    print(f'Geopotential:   U = {U}')
    print(f'R (Sun):    {sun.disturbing_function(satellite) / U}')
    print(f'R (Moon):   {moon.disturbing_function(satellite) / U}')


def triangulation():
    """Определение положения спутника по наблюдениям из двух пунктов методом триангуляции"""
    
    point1 = Point(
        (-3811.723938399338, -2775.741899662411, -3286.0884872272504))
    point2 = Point((225.87627253590838, 5477.069490937063, 2938.8013951824414))

    gamma1 = Angles(decimal=36.046824810480096)
    gamma2 = Angles(decimal=267.635619629101)

    delta1 = Angles(decimal=34.908846560020685)
    delta2 = Angles(decimal=-28.159853782057976)

    triang = Triangulate.Triangulate(point1, point2,
                                     gamma1, delta1,
                                     gamma2, delta2)
    (x1, y1, z1), (x2, y2, z2) = triang.get_coords()

    print(f'x\'={x1}        x"={x2}')
    print(f'y\'={y1}        y"={y2}')
    print(f'z\'={z1}        z"={z2}')


def get_recieve_coords():
    coords = np.array([
        [-8420.491146597293, -5468.7108279490885, -3225.2768744335112],
        [10627.866333593165, 12197.779078296462, 12749.899708981226],
        [3028.8049577894853, -7254.328380486594, 10674.147072829894],
        [-3394.545889436539, 12134.310645747004, -8333.83749367311],
        [-15647.953138265704, -15866.541879902587, -18663.162790105933],
        [11884.641571123015, 2391.714901148104, 5756.097816655203],
        [-6126.52417257099, -14125.496799667424, 1975.0853234620217],
        [19834.743174580144, -16057.9045301084, -3613.1923164238387],
        [-4751.389861824661, 6356.270245925387, -1419.8451589415224],
        [-7781.245391093064, -5945.1160422309795, -16763.688018102195],
        [1040.1788489550577, -18635.119624119056, -17415.12505575222],
        [-11664.034164720966, 9291.433375198221, -11833.343574028651],
    ])

    d = np.array([
        5216.459164383287,
        27878.259759555014,
        19138.890127546383,
        14041.489343059386,
        22112.067298510538,
        21221.54259758233,
        14671.044163276372,
        29740.46543313604,
        8601.37945433659,
        12641.922223583239,
        22351.766724051802,
        13771.081214445134
    ])

    dt = 0

    q0 = 1 / len(d) * sum(coords)
    q0 = np.append(q0, dt)

    dif = 1

    while np.linalg.norm(dif) > 1e-6:
        dt = q0[3]

        rho = Triangulate.get_rho(XYZ=q0, xyz=coords)
        delta = rho - d - Mechanics.c * dt

        F = Triangulate.F(XYZ=q0, xyz=coords.T, rho=rho, delta=delta)
        dF = Triangulate.dF(XYZ=q0, xyz=coords.T, rho=rho, delta=delta)

        q0 -= np.linalg.inv(dF) @ F
        dif = F

    print(
        f'X = {float(q0[0])}\nY = {float(q0[1])}\nZ = {float(q0[2])}\ndt = {float(q0[3])}')

    return


def satellite_leveling():
    """Спутниковое нивелирование"""
    
    # Вариант:           8
    # Координаты спутника, км:
    # 12419.4473703004422553        -16144.7191515973447915        -1879.43864706130342057
    # Время распространения сигнала, мс:   93.8845680925201899678
    xs = 12419.4473703004422553
    ys = -16144.7191515973447915
    zs = -1879.43864706130342057
    dt = 93.8845680925201899678 / 1000

    leveling = Leveling.Leveling(xs=xs, ys=ys, zs=zs)

    lmd = leveling.Newton_method()

    x, y, z = leveling.coords(lmd=lmd)

    difs = [x-xs, y-ys, z-zs]
    d = Math.radius(difs)

    h = d - Mechanics.c * dt/2

    print(f'h = {h} km')


if __name__ == "__main__":
    # disturb()
    # triangulation()
    # get_recieve_coords()
    satellite_leveling()
