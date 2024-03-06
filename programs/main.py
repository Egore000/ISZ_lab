import numpy as np

from config import *
from tools import Filer, Grapher, Parser
from mechanics import Mechanics
from MathPy import Math
from objects import Satellite, Earth, Sun, Moon, Point


def orbit(type: str):
    earth = Earth(EARTH_PATH)
    satellite = Satellite(type=type, t0=0)

    evolution = np.array(list(satellite.evolution())).T
    coords = np.array(list(satellite.coords)).T
    
    grapher = Grapher(custom_rcParams, projection='3d')
    grapher.print(evolution, c='black')
    grapher.print(earth.coords.T)

    grapher.print(coords, c='red', s=10)
    grapher.fig.suptitle(type + ' спутник')
    grapher.show()
    

def route(type: str):
    satellite = Satellite(type=type)
    route = satellite.route(date)
    
    grapher = Grapher(custom_rcParams, projection=None)
    grapher.print(route, title=type + ' спутник', **marker)
    grapher.show()


def animation(type: str):
    satellite = Satellite(type=type)
    earth = Earth(EARTH_PATH)
    grapher = Grapher(projection='3d')
    grapher.animation(satellite, earth, save=0, title='Резонанс')


def disturb():
    sun = Sun((
             -1.848432950986775E-01 * Mechanics.ae, 
              8.860598024929441E-01 * Mechanics.ae, 
              3.841027854355354E-01 * Mechanics.ae
              ))

    moon = Moon((
                (-1.863831425711169E-01 + 1.848432950986775E-01) * Mechanics.ae,
                ( 8.878791230277183E-01 - 8.860598024929441E-01) * Mechanics.ae,
                ( 3.850801888852833E-01 - 3.841027854355354E-01) * Mechanics.ae
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
    point1 = Point((-3811.723938399338, -2775.741899662411, -3286.0884872272504))
    point2 = Point((225.87627253590838, 5477.069490937063, 2938.8013951824414))
    
    gamma1 = Angles(decimal=36.046824810480096) 
    gamma2 = Angles(decimal=267.635619629101)

    delta1 = Angles(decimal=34.908846560020685)
    delta2 = Angles(decimal=-28.159853782057976)

    (x1, y1, z1), (x2, y2, z2) = Math.Triangulate.get_coords(point1, point2,
                                                            gamma1, delta1, 
                                                            gamma2, delta2)
    print(f'x\'={x1}        x"={x2}')
    print(f'y\'={y1}        y"={y2}')
    print(f'z\'={z1}        z"={z2}')


if __name__=="__main__":
    # route('Тестовый')
    # orbit('Тестовый')
    # animation('Тестовый')
    # route('Геостационарный')
    # orbit('Геостационарный')
    # animation('Геостационарный')

    # disturb()

    triangulation()