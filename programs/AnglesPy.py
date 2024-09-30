import math
from typing import Any

k = 180/math.pi

class Angles:
    '''
    Класс, прeдназначенный для облегчения работы с тригонометрическими углами.
    
    Для начала необходимо инициализировать угол, передав в качестве аргументов
    нужные значения угловых градусов, минут и секунд, либо указать значение 
    в радианах или десятичных градусах:
    >>> angle = AnglesPy.Angles(degree=30, minute=25, second=10)
    # <Angles: 30° 25' 10">
    >>> angle = AnglesPy.Angles(17, 35, 11)
    # <Angles: 17° 35' 11">
    >>> angle = AnglesPy.Angles(rad=math.pi)
    # <Angles: 180° 00' 00">
    >>> angle = AnglesPy.Angles(decimal=45.5)
    # <Angles: 45° 30' 00">

    Параметры:
    ----------
    `deg` : `int` - Градусы
    >>> angle.deg
    30

    `min` : `int` - Минуты
    >>> angle.min
    25

    `sec` : `float` - Секунды
    >>> angle.sec
    10

    `decimal` : `float` - Угол в десятичных градусах
    >>> angle.decimal
    30.419444444444444

    `rad` : `float` - Радианная мера угла
    >>> angle.rad
    0.5309194621830529

    `dms` : `str` - Строка в формате d° m' s" для удобного вывода на консоль.
    >>> angle.dms
    "30° 25' 10.00\""

    Для ввода угла в радианах или десятичных градусах используются
    параметры `rad` или `decimal`: 
    >>> angle = AnglesPy.Angles(rad=3.14)
    >>> angle = AnglesPy.Angles(decimal=14.5)
    '''
    params = {
    'units': 'rad',
    'precision': 2,
    }

    def __init__(self, degree=0, minute=0, second=0, rad=None, decimal=None):
        if rad or decimal:
            if rad:
                self.rad = rad
                self.decimal = rad * k
            elif decimal:
                self.decimal = decimal
                self.rad = decimal / k
            self.deg, self.min, self.sec = self.__get_dms(self.decimal)
        else:
            if degree < 0:
                minute = -abs(minute)
                second = -abs(second)
            elif degree == 0 and minute < 0:
                second = -abs(second)

            if second >= 60:
                minute += second // 60
                second %= 60 

            if minute >= 60:
                degree += minute // 60
                minute %= 60

            self.deg = degree
            self.min = minute
            self.sec = second

            self.decimal = self.deg + self.min/60 + self.sec/3600
            self.rad = self.decimal / k

        precision = self.params.get('precision')
        sec = abs(self.sec)

        self.dms = f"{self.deg}° {abs(self.min)}' {sec:.{precision}f}\""
        if self.deg == 0:
            if self.min == 0:
                if self.sec < 0:
                    self.dms = '-' + self.dms
            elif self.min < 0:
                self.dms = '-' + self.dms

    def __add__(self, other):
        if isinstance(other, Angles):
            return Angles(rad=self.rad + other.rad)
        else:
            return Angles(rad=self.rad + other)
       
    def __sub__(self, other):
        if isinstance(other, Angles):
            return Angles(rad=self.rad - other.rad)
        else:
            return Angles(rad=self.rad - other)

    def __mul__(self, other):
        if isinstance(other, Angles):
            return self.rad * other.rad
        else:
            return Angles(rad=self.rad * other)
    
    def __truediv__(self, other):
        if isinstance(other, Angles):
            return self.rad / other.rad
        else:
            return Angles(rad=self.rad / other)
        
    def __radd__(self, other):
        return self + other
    
    def __rsub__(self, other):
        return self - other
    
    def __rmul__(self, other):
        return self * other

    def __iadd__(self, other):
        return self + other

    def __isub__(self, other):
        return self - other

    def __imul__(self, other):
        return self * other

    def __itruediv__(self, other):
        return self / other

    def __eq__(self, other):
        if isinstance(other, Angles):
            return self.rad == other.rad
        else:
            return self.rad == other

    def __ne__(self, other):
        if isinstance(other, Angles):
            return self.rad != other.rad
        else:
            return self.rad != other

    def __lt__(self, other):
        if isinstance(other, Angles):
            return self.rad < other.rad
        else:
            return self.rad < other

    def __gt__(self, other):
        if isinstance(other, Angles):
            return self.rad > other.rad
        else:
            return self.rad > other

    def __le__(self, other):
        if isinstance(other, Angles):
            return self.rad <= other.rad
        else:
            return self.rad <= other

    def __ge__(self, other):
        if isinstance(other, Angles):
            return self.rad >= other.rad
        else:
            return self.rad >= other

    def __float__(self):
        return self.rad

    def __str__(self):
        return self.dms

    def __abs__(self) -> float:
        return abs(self.rad)
    
    def __neg__(self):
        return Angles(rad=-self.rad)
    
    def __repr__(self):
        precision = self.params.get("precision")
        return repr(f'<Angles: {self.deg}° {self.min}` {self.sec:.{precision}f}">')
    
    def reduce(self):
        value = self.rad - self.rad // math.pi * math.pi
        return Angles(rad=value)
    
    @staticmethod
    def __get_dms(angle: Any) -> tuple[int, int, Any]:
        deg = int(angle)
        minute = int((angle - deg) * 60)
        sec = (angle - deg - minute/60) * 3600
        return (deg, minute, sec)


def main():
    return
    
if __name__=="__main__":
    main()