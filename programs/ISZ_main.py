import numpy as np

from config import *
from MathPy import Math
from tools import Filer, Grapher, Parser
from objects import Satellite, Observer, Earth, GLONASS
from mechanics import Mechanics


__doc__ = """
Задачи по исследованию динамики ИСЗ
"""

def orbit(type: str):
    """Построение орбиты спутника в 3D"""

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
    """Построение трассы спутника"""

    satellite = Satellite(type=type)
    route = satellite.route(date)
    initial_coords = route[0]

    grapher = Grapher(custom_rcParams, projection=None)
    grapher.print(route, title=type + ' спутник', **marker)
    grapher.print(initial_coords, c='red')
    grapher.show()


def animation(type: str):
    """Анимированный полёт спутника вокруг Земли в 3D"""

    satellite = Satellite(type=type)
    earth = Earth(EARTH_PATH)
    grapher = Grapher(projection='3d')
    grapher.animation(satellite, earth, save=1, title='Резонанс')


def glonass(time):
    """Построение трасс и орбит спутников ГЛОНАСС"""
    g = GLONASS()

    # g.current_position()
    g.get_positions(time)
    g.get_routes(time)
    g.get_orbits(time)


def visibility_range(time):
    """Вычисление зон видимости спутников ГЛОНАСС"""
    glonass = GLONASS()

    jd = Math.get_JD(time)
    t = jd * 86400

    geo = Satellite(type='Геостационарный')
    geo.coords, _ = Mechanics.get_coords(geo, t)

    satellites = glonass.get_satellites(time)

    in_zone = []
    for sat in satellites:
        R = Math.get_perpendicular(sat, geo)

        print(sat.type, R)

        if R > Earth.Radius:
            in_zone.append(sat.type)

    print(in_zone)
    glonass.visibility_areas(time)
    # glonass.get_orbits(time)


def is_satellite_in_visibility_range_for_observer(satellite: Satellite, observer: Observer, time: str) -> bool:
    """Проверка видимости спутника из точки наблюдения"""

    H = Math.get_sidereal_time(time)
    LST = observer.local_sidereal_time(sidereal_time=H)
    observer_x, observer_y, observer_z = observer.coords_CRS(sidereal_time=H)

    r_xyz = np.array([
        satellite.x - observer_x,
        satellite.y - observer_y,
        satellite.z - observer_z,
    ])

    r_sez = Mechanics.HTS(r_xyz, local_sidereal_time=LST, lat=observer._lat)
    r_s, r_e, r_z = r_sez

    r = Math.radius(r_sez)

    h = np.arcsin(r_z / r)
    return h > np.pi/6


def visibility_areas(time: str):
    """Построение зон видимости"""

    def init_observer(lat, lon, height, time):
        observer = Observer(lat=lat, lon=lon, height=height)
        H = Math.get_sidereal_time(time)
        LST = observer.local_sidereal_time(sidereal_time=H)
        return observer, LST

    lat0 = Angles()
    lon0 = Angles()

    lat_Tomsk = Angles(56, 29, 19)
    lon_Tomsk = Angles(84, 57, 8)

    height = 0

    observer0, LST0 = init_observer(lat0, lon0, height, time)
    observer_Tomsk, LST = init_observer(lat_Tomsk, lon_Tomsk, height, time)

    h_min = Angles(30)
    Azimuth = [Angles(deg) for deg in range(1, 361)]
    a = 23047.00953687991
    r = (Earth.Radius + np.sqrt(4 * a**2 - 3 * Earth.Radius**2)) / 2

    grapher = Grapher(custom_rcParams=custom_rcParams, projection=None)
    grapher.print(observer0.coords, c='green', s=100)
    grapher.print(observer_Tomsk.coords, c='red', s=100)

    for Az in Azimuth:
        r_s = r * Math.cos(Az) * Math.cos(h_min)
        r_e = -r * Math.sin(Az) * Math.cos(h_min)
        r_z = r * Math.sin(h_min)

        r_sez = np.array([r_s, r_e, r_z])

        r0 = Mechanics.from_HTS_to_CRS(
            vector=r_sez, local_sidereal_time=LST0, lat=observer0._lat)
        r_Tomsk = Mechanics.from_HTS_to_CRS(
            vector=r_sez, local_sidereal_time=LST, lat=observer_Tomsk._lat)

        lmd0, phi0 = Mechanics.from_CRS_to_TRS(vector=r0, time=time)
        lmd_Tomsk, phi_Tomsk = Mechanics.from_CRS_to_TRS(
            vector=r_Tomsk, time=time)

        grapher.print(data=(lmd0, phi0), c='green', s=1)
        grapher.print(data=(lmd_Tomsk, phi_Tomsk), c='red', s=1)

    glonass = GLONASS()
    sats = glonass.get_satellites(time)

    vis0 = []
    vis_Tomsk = []
    for sat in sats:
        coords = Math.get_lmd_phi(sat.coords)
        grapher.print(coords, c='k', s=20)
        grapher.ax.annotate(int(sat.type), (coords[0].decimal,
                                            coords[1].decimal + 2))

        if is_satellite_in_visibility_range_for_observer(sat, observer0, time):
            vis0.append(sat.type)

        if is_satellite_in_visibility_range_for_observer(sat, observer_Tomsk, time):
            vis_Tomsk.append(sat.type)

    print(f'0 = {vis0}', f'Tomsk = {vis_Tomsk}', sep='\n')

    grapher.ax.set_ylim(-90, 90)
    grapher.ax.set_xlim(-180, 180)
    grapher.show()
    return


def trisection():
    """Определение координат спутника методом трисекции"""

    coords = np.array([
        [  48853.4796304855, -1157502.1971707679, -229687.1713272969],
        [  32862.01144554671,  157502.1971707679, -208608.9455096204],
        [-167247.2391357287,  -158824.7863056938, -127309.0503967787]
    ])
    satellite = Satellite(type='Геостационарный')

    rho = np.array([
        Math.radius(np.array(satellite.coords) - coords[:: ,0]),
        Math.radius(np.array(satellite.coords) - coords[::, 1]),
        Math.radius(np.array(satellite.coords) - coords[::, 2]),
    ])

    x, y, z = 0, 0, 0

    eps = 1e-9
    dif = 1
    iteration = 0
    while dif >= eps:
        iteration += 1

        _x = coords[0, ::]
        _y = coords[1, ::]
        _z = coords[2, ::]

        A = np.array([
            [
                (x - _x[0]) / rho[0],  
                (y - _y[0]) / rho[0], 
                (z - _z[0]) / rho[0]
            ],
            [
                (x - _x[1]) / rho[1],
                (y - _y[1]) / rho[1],
                (z - _z[1]) / rho[1]
            ],
            [
                (x - _x[2]) / rho[2],
                (y - _y[2]) / rho[2],
                (z - _z[2]) / rho[2]
            ]
        ])

        xyz = np.array([x, y, z])
        _rho = np.array([
            Math.radius(xyz - coords[::, 0]),
            Math.radius(xyz - coords[::, 1]),
            Math.radius(xyz - coords[::, 2]),
        ])

        delta_rho = rho - _rho

        delta_x = np.linalg.inv(A) @ delta_rho
        xyz = xyz + delta_x
        x, y, z = xyz

        print(f'{iteration}: {x=} км    {y=} км    {z=} км')

        dif = Math.radius(delta_x)
        
    print(f'Координаты спутника: {satellite.coords}',
        f'Полученные координаты: {x, y, z}', 
        f'Количество итераций: {iteration}', 
        sep='\n')
    


if __name__ == "__main__":
    # route('Тестовый')
    # orbit('Тестовый')
    # animation('Тестовый')
    # route('Геостационарный')
    # orbit('Геостационарный')
    animation('Геостационарный')
    # glonass(date)
    # visibility_range(date)
    # visibility_areas(date)
    # trisection()