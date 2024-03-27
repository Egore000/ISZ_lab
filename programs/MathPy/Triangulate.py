import numpy as np
import numpy.typing as npt
from math import cos, sin

from AnglesPy import Angles
from . import Math
from mechanics import Mechanics

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


def _Fx(X: float, 
        x: list[float], 
        rho: list[float],
        delta: list[float]) -> float:
    
    '''F_x'''
    sum = 0
    for i in range(len(x)):
        sum += delta[i] * (X - x[i]) / rho[i]
    return 2 * sum

def _Ft(delta: list[float]) -> float:
    '''F_t'''
    return -2 * Mechanics.c * sum(delta)

def dF_diagonal(X: float,
         x: list[float],
         rho: list[float],
         delta: list[float]) -> float:
    
    sum = 0
    for i in range(len(x)):
        sum += ( (X - x[i])**2 + delta[i] * (rho[i] - (X - x[i])**2 / rho[i]) ) / rho[i]**2

    return 2 * sum

def dF_xy(XY: list[float],
          xy: list[list[float]],
          rho: list[float],
          delta: list[float]) -> float:
    
    X, Y = XY
    x, y = xy

    sum = 0
    for i in range(len(x)):
        sum += (X - x[i]) * (Y - y[i]) * (rho[i] - delta[i]) / rho[i]**3
    return 2 * sum


def dF_xt(X: float,
          x: list[float],
          rho: list[float]) -> float:

    sum = 0
    for i in range(len(x)):
        sum += (X - x[i]) / rho[i]
    
    return -2 * Mechanics.c * sum


def F(XYZ: list[float],
      xyz: list[list[float]],
      rho: list[float],
      delta: list[float]) -> npt.NDArray:
    
    X, Y, Z, _ = XYZ
    x, y, z = xyz

    lst = [
        _Fx(X=X, x=x, rho=rho, delta=delta),
        _Fx(X=Y, x=y, rho=rho, delta=delta),
        _Fx(X=Z, x=z, rho=rho, delta=delta),
        _Ft(delta=delta)
    ]
    return np.array(lst) 


def dF(XYZ: list[float],
       xyz: list[list[float]],
       rho: list[float],
       delta: list[float]) -> npt.NDArray:
    
    X, Y, Z, _ = XYZ
    x, y, z = xyz

    params = {
        'rho': rho,
        'delta': delta         
    }

    matrix = np.array([
        [
            dF_diagonal(X=X, x=x, **params),  
            dF_xy(XY=(X, Y), xy=[x, y], **params), 
            dF_xy(XY=(X, Z), xy=[x, z], **params),
            dF_xt(X=X, x=x, rho=rho)   
        ],
        [
            dF_xy(XY=(X, Y), xy=[x, y], **params), 
            dF_diagonal(X=Y, x=y, **params),  
            dF_xy(XY=(Y, Z), xy=[y, z], **params),
            dF_xt(X=Y, x=y, rho=rho)   
        ],
        [
            dF_xy(XY=(X, Z), xy=[x, z], **params), 
            dF_xy(XY=(Y, Z), xy=[y, z], **params),
            dF_diagonal(X=Z, x=z, **params),  
            dF_xt(X=Z, x=z, rho=rho)   
        ],   
        [
            dF_xt(X=X, x=x, rho=rho),
            dF_xt(X=Y, x=y, rho=rho),  
            dF_xt(X=Z, x=z, rho=rho),  
            2 * len(x) * Mechanics.c**2
        ]
    ])
    return matrix

def get_rho(XYZ: list[float],
            xyz: list[list[float]]) -> list[float]:
    
    X, Y, Z, _ = XYZ
    lst = []
    for item in xyz:
        x, y, z = item
        lst.append(Math.radius(((X - x), (Y - y), (Z - z))))
    return lst