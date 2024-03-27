import numpy as np

from config import *
from tools import Filer, Grapher, Parser
from objects import Satellite, Earth, GLONASS


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

def glonass():
    g = GLONASS()
    # g.current_position()
    g.get_positions('27.03.2024 08:19:40')
    # g.get_routes('14.03.2024 11:58:04')
    # g.get_orbits('14.03.2024 10:16:34')


if __name__=="__main__":
    # route('Тестовый')
    # orbit('Тестовый')
    # animation('Тестовый')
    # route('Геостационарный')
    # orbit('Геостационарный')
    # animation('Геостационарный')
    glonass()