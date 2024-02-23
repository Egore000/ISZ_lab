import numpy as np

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