import numpy as np
from math import cos, sin

from AnglesPy import Angles


class Math:
    '''
    Класс с основными математическими вычислениями
    для динамики ИСЗ
    '''
    @staticmethod
    def factorial(x: int) -> int:
        '''Вычисление факториала целого числа `x`'''
        if not isinstance(x, int):
            raise ValueError('n!: n must be int')
        result = 1
        for i in range(1, x):
            result *= i
        return result
    
    @staticmethod
    def LPn(n: int, z: float) -> float:
        '''
        Функция Лежандра `Pn(z)`
        ```
        -1 <= z <= 1
        ```
        '''
        if not -1 <= z <= 1:
            raise ValueError('z must be in [-1, 1]')
        if n == 0:
            return 1
        if n == 1:
            return z
        return 1/n * ((2*n - 1) * z * Math.LPn(n-1, z) - (n - 1) * Math.LPn(n-2, z))
    
    @staticmethod
    def LPnm(n: int, m: int, z: float) -> float:
        '''
        Присоединённая функция Лежандра `Pnm(z)`
        ```
        -1 <= z <= 1
        ```
        '''
        if not -1 <= z <= 1:
            raise ValueError('z is not between -1 and 1')
        if m == 0 and n < 2:
            if n == 0: 
                return 1
            elif n == 1: 
                return z
        elif m == n:
            return (2*n - 1) *  np.sqrt(1 - z**2) * Math.LPnm(n-1, n-1, z)
        else:
            if n - 1 < m:
                return 0
            elif n - 2 < m:
                return (2*n - 1) * z * Math.LPnm(n-1, m, z) / (n-m)
            else:
                return ((2*n - 1) * z * Math.LPnm(n-1, m, z) - (n - 1 + m) * Math.LPnm(n-2, m, z)) / (n-m)

    @staticmethod
    def Norm_LPnm(n: int, m: int, z: float) -> float:
        '''
        Нормированная функция Лежандра 
        ```
        -1 <= z <= 1
        ```
        '''
        coef = np.sqrt(2*n + 1) * np.sqrt(2 * Math.factorial(n - m) / Math.factorial(n + m))
        return coef * Math.LPnm(n, m, z)
    
    @staticmethod
    def get_lmd_phi(args: tuple[float]) -> tuple[float]:
        '''
        Вычисление географических координат по декартовым
        '''
        (x, y, z) = args
        lmd = Angles(rad=np.arctan2(y, x))
        phi = Angles(rad=np.arctan(z / np.sqrt(x**2 + y**2)))
        return (lmd, phi)
        
    @staticmethod
    def radius(args: tuple[float]) -> float:
        '''Радиус-вектор'''
        (x, y, z) = args
        return np.sqrt(x**2 + y**2 + z**2)

    @staticmethod
    def sid2000(jd: float) -> float:
        '''
        Вычисление звёздного времени sid2000
            `jd` - юлианская дата
        '''
        jd2000 = 2451545
        jdyear = 36525
        m = jd - int(jd) - 0.5
        d = jd - m - jd2000
        t = (d + m)/jdyear
        mm = m * 86400
        s = (24110.54841 + mm + 236.555367908 * (d + m) + (0.093104 * t - 6.21E-6 * t**2) * t) / 86400 * 2*np.pi
        return s
    
    @staticmethod
    def get_JD(date: str) -> float:
        '''
        Вычисление юлианской даты
        '''
        year, month, day = Math.get_daytime(date)
        
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
    
    @staticmethod
    def get_daytime(date: str) -> tuple[int, int, float]:
        '''
        Приведение даты в строков формате к виду `(year, month, day)`
            `year`: `int`
            `month`: `int`
            `day`: `float`
        '''
        date, time = date.split()
        day, month, year = map(int, date.split('.'))
        hour, min = map(int, time.split(':'))
        day += (hour + min/60)/24
        return (year, month, day)
    
    class Triangulate:
        '''
        Подкласс для работы с методами триангуляции
        '''
        def __init__(self, point1, point2, 
                    gamma1: Angles, delta1: Angles,
                    gamma2: Angles, delta2: Angles):
            self.point1 = point1
            self.point2 = point2
            self.g1 = gamma1
            self.g2 = gamma2
            self.d1 = delta1
            self.d2 = delta2

        @staticmethod
        def get_components(gamma: Angles, delta: Angles) -> tuple[float]:
            l = cos(gamma) * cos(delta)
            m = sin(gamma) * cos(delta)
            n = sin(delta)
            return (l, m, n)

        def get_coords(self) -> tuple[float]:
            x1, y1, z1 = self.point1.coords
            x2, y2, z2 = self.point2.coords

            dx, dy, dz = x2-x1, y2-y1, z2-z1

            l1, m1, n1 = Math.Triangulate.get_components(self.g1, self.d1)
            l2, m2, n2 = Math.Triangulate.get_components(self.g2, self.d2)

            cos_a = l1 * l2 + m1 * m2 + n1 * n2
            sin_a = np.sqrt(1 - cos_a ** 2)

            F1 = dx*l1 + dy*m1 + dz*n1
            F2 = dx*l2 + dy*m2 + dz*n2

            rho1 = (F1 - F2*cos_a) / sin_a**2
            rho2 = (F1*cos_a - F2) / sin_a**2

            x1 = self.point1.x + rho1 * l1
            x2 = self.point2.x + rho2 * l2

            y1 = self.point1.y + rho1 * m1
            y2 = self.point2.y + rho2 * m2

            z1 = self.point1.z + rho1 * n1
            z2 = self.point2.z + rho2 * n2

            return (x1, y1, z1), (x2, y2, z2)