import math as m

from AnglesPy import Angles

class Ellipsoid:
    def __init__(self, a: float, b: float):
        self.ecc_2 = 1 - (b / a)**2
        self.ecc2_2 = (a / b)**2 - 1
        self.ecc = m.sqrt(self.ecc_2)
        self.ecc2 = m.sqrt(self.ecc2_2)

    def N(self, B: Angles) -> float:
        return self.a * (1 - self.ecc_2 * m.sin(B)**2)**(-0.5)
    
    def M(self, B: Angles) -> float:
        return self.a * (1 - self.ecc_2) * (1 - self.ecc_2 * m.sin(B)**2)**(-1.5)
    
    def R(self, B: Angles) -> float:
        return m.sqrt(self.M(B) * self.N(B))
    

class Krasovsky(Ellipsoid):
    def __init__(self):
        self.a = 6378.245
        self.b = 6356.863019
        self.f = 1 / 298.3
        super().__init__(self.a, self.b)


class WGS84(Ellipsoid):
    def __init__(self):
        self.a = 6378.137
        self.b = 6356.752
        self.f = 1 / 298.257
        super().__init__(self.a, self.b)


krasovsky = Krasovsky()
wgs84 = WGS84()