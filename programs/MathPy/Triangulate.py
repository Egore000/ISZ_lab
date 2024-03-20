import numpy as np
from math import cos, sin

from AnglesPy import Angles
from . import Math


class Triangulate:
    '''
    Класс для работы с методами триангуляции
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

        l1, m1, n1 = Triangulate.get_components(self.g1, self.d1)
        l2, m2, n2 = Triangulate.get_components(self.g2, self.d2)

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