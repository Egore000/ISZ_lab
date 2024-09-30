from GeoPY import krasovsky

a = krasovsky.a
b = krasovsky.b


class Leveling:
    """
    Нивелирование
    """

    def __init__(self, xs, ys, zs):
        self.xs = xs
        self.ys = ys
        self.zs = zs

    def coords(self, lmd: float) -> tuple[float]:
        """Получение координат на эллипсоиде"""

        x = a**2 / (a**2 + lmd) * self.xs
        y = a**2 / (a**2 + lmd) * self.ys
        z = b**2 / (b**2 + lmd) * self.zs
        return x, y, z

    def derivative_coords(self, lmd: float) -> tuple[float]:
        """Производные от координат"""
        
        dx = - (a / (a**2 + lmd))**2 * self.xs
        dy = - (a / (a**2 + lmd))**2 * self.ys
        dz = - (b / (b**2 + lmd))**2 * self.zs
        return dx, dy, dz

    @staticmethod
    def _f(x: float, y: float, z: float):
        """f(lmd). Уравнение эллипсоида"""
        return (x/a)**2 + (y/a)**2 + (z/b)**2 - 1

    @staticmethod
    def _df(x, y, z,
            dx, dy, dz):
        """f'(lmd)"""
        return 2 * (x*dx / a**2 + y*dy / a**2 + z*dz / b**2)

    def Newton_method(self, eps=1e-6) -> float:
        """Метод Ньютона для вычисления долготы"""

        lmd = 0
        dif = 1
        while abs(dif) > eps:
            x, y, z = self.coords(lmd=lmd)
            dx, dy, dz = self.derivative_coords(lmd=lmd)

            dif = self._f(x, y, z) / self._df(x, y, z, dx, dy, dz)
            lmd -= dif

        return lmd
