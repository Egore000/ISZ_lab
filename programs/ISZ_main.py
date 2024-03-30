import numpy as np

from config import *
from MathPy import Math
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
    initial_coords = route[0]

    grapher = Grapher(custom_rcParams, projection=None)
    grapher.print(route, title=type + ' спутник', **marker)
    grapher.print(initial_coords, c='red')
    grapher.show()


def animation(type: str):
    satellite = Satellite(type=type)
    earth = Earth(EARTH_PATH)
    grapher = Grapher(projection='3d')
    grapher.animation(satellite, earth, save=0, title='Резонанс')


def glonass(time):
    g = GLONASS()

    # g.current_position()
    g.get_positions(time)
    g.get_routes(time)
    g.get_orbits(time)


if __name__ == "__main__":
    # route('Тестовый')
    # orbit('Тестовый')
    # animation('Тестовый')
    # route('Геостационарный')
    # orbit('Геостационарный')
    # animation('Геостационарный')
    glonass('14.03.2024 11:58:04')
