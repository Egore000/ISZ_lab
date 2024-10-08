import numpy as np
from math import cos, sin

from AnglesPy import Angles


def factorial(x: int) -> int:
    """Вычисление факториала целого числа `x`"""
    if not isinstance(x, int):
        raise ValueError("n!: n must be int")
    result = 1
    for i in range(1, x):
        result *= i
    return result


def LPn(n: int, z: float) -> float:
    """
    Функция Лежандра `Pn(z)`
    ```
    -1 <= z <= 1
    ```
    """
    if not -1 <= z <= 1:
        raise ValueError("z must be in [-1, 1]")
    if n == 0:
        return 1
    if n == 1:
        return z
    return 1/n * ((2*n - 1) * z * LPn(n-1, z) - (n - 1) * LPn(n-2, z))


def LPnm(n: int, m: int, z: float) -> float:
    """
    Присоединённая функция Лежандра `Pnm(z)`
    ```
    -1 <= z <= 1
    ```
    """
    if not -1 <= z <= 1:
        raise ValueError("z is not between -1 and 1")
    if m == 0 and n < 2:
        if n == 0:
            return 1
        elif n == 1:
            return z
    elif m == n:
        return (2*n - 1) * np.sqrt(1 - z**2) * LPnm(n-1, n-1, z)
    else:
        if n - 1 < m:
            return 0
        elif n - 2 < m:
            return (2*n - 1) * z * LPnm(n-1, m, z) / (n-m)
        else:
            return ((2*n - 1) * z * LPnm(n-1, m, z) - (n - 1 + m) * LPnm(n-2, m, z)) / (n-m)


def Norm_LPnm(n: int, m: int, z: float) -> float:
    """
    Нормированная функция Лежандра 
    ```
    -1 <= z <= 1
    ```
    """
    coef = np.sqrt(2*n + 1) * np.sqrt(2 * factorial(n - m) / factorial(n + m))
    return coef * LPnm(n, m, z)


def get_lmd_phi(args: tuple[float]) -> tuple[Angles, Angles]:
    """
    Вычисление географических координат по декартовым
    """
    (x, y, z) = args
    lmd = Angles(rad=np.arctan2(y, x))
    phi = Angles(rad=np.arctan(z / np.sqrt(x**2 + y**2)))
    return (lmd, phi)


def radius(args: tuple[float]) -> float:
    """Радиус-вектор"""
    (x, y, z) = args
    return np.sqrt(x**2 + y**2 + z**2)


def sid2000(jd: float) -> float:
    """
    Вычисление звёздного времени sid2000
        `jd` - юлианская дата
    """
    jd2000 = 2451545
    jdyear = 36525
    m = jd - int(jd) - 0.5
    d = jd - m - jd2000
    t = (d + m)/jdyear
    mm = m * 86400
    s = (24110.54841 + mm + 236.555367908 * (d + m) +
         (0.093104 * t - 6.21E-6 * t**2) * t) / 86400 * 2*np.pi
    return s


def get_JD(date: str) -> float:
    """Вычисление юлианской даты"""
    
    year, month, day = get_daytime(date)

    date = year + month/100 + day/1e4
    i = int(date)
    m1 = (date - i) * 100
    me = int(m1)
    d = (m1 - me) * 100
    if me <= 2:
        i -= 1
        me += 12
    jd = int(365.25 * i) + int(30.6001 * (me + 1)) + d + 1720994.5

    if date < 1582.1015:
        return jd
    ja = int(i/100)
    jb = 2 - ja + int(ja/4)
    return jd + jb


def get_daytime(date: str) -> tuple[int, int, float]:
    """
    Приведение даты в строков формате к виду `(year, month, day)`
        `year`: `int`
        `month`: `int`
        `day`: `float`
    """
    date, time = date.split()
    day, month, year = map(int, date.split("."))
    hour, min, sec = map(int, time.split(":"))
    day += (hour + min/60 + sec/3600)/24
    return (year, month, day)


def get_perpendicular(satellite1, satellite2) -> float:
    x12 = satellite1.coords[0] - satellite2.coords[0]
    y12 = satellite1.coords[1] - satellite2.coords[1]
    z12 = satellite1.coords[2] - satellite2.coords[2]

    r12 = radius((x12, y12, z12))

    r1 = radius(satellite1.coords)
    r2 = radius(satellite2.coords)

    return np.sqrt((2 * r12 * r1)**2 - (r12**2 + r1**2 - r2**2)**2) / (2 * r12)


def get_sidereal_time(date: str) -> Angles:
    """Звёздное время"""
    _, _, day = get_daytime(date)

    sidereal_time = (day - int(day)) / 24 * 360
    
    return Angles(decimal=sidereal_time)