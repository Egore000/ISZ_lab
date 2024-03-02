import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from config import * 
from MathPy import *


class Object:
    def __init__(self, coords: tuple):
        self.coords = coords

    def disturbing_function(self, satellite):
        (x, y, z) = satellite.coords
        (x_, y_, z_) = self.coords

        dx = (x-x_, y-y_, z-z_)
        r_ = Math.radius(self.coords)
        delta = Math.radius(dx)

        return self.fm * (1 / delta - (x*x_ + y*y_ + z*z_) / r_**3)


class Earth:
    Mass = 2e24
    Radius = 6378
    fm = 398600.5
    def __init__(self, file: str):
        self.coords = Filer.read(file)


class Sun(Object):
    fm = 332946.04877304828 * Earth.fm
    

class Moon(Object):
    fm = Earth.fm / 81.30056822149722
    


class Filer:
    '''Класс для чтения и записи в файлы'''
    @staticmethod
    def read(file: str) -> list:
        '''Чтение данных из файла `file`'''
        data = []
        with open(file, 'r') as f:
            for line in f:
                line_filtered = [float(x) for x in line.strip().split(' ') if x != '']
                data.append(np.array(line_filtered))
        return np.array(data)

    @staticmethod
    def read_coefficients(path):
        '''Чтение данных с коэффициентами геопотенциала'''
        array = {}
        with open(path, 'r') as file:
            for i, line in enumerate(file):
                if not i:
                    continue
                n, m, Cnm, Snm = line.split()
                n = int(n)
                m = int(m)
                Cnm = float(Cnm)
                Snm = float(Snm)
                
                array[(n, m)] = (Cnm, Snm)
        return array
    
    @staticmethod
    def read_file(path) -> tuple[list[float] | float]:
        '''
        Чтение данных из path

        Возрвращает следующие данные::
        -----------
            `time` - время

            `coords` - координаты
            
            `velocities` - скорости
            
            `megno` - параметр MEGNO
            
            `mean_megno` - осреднённый MEGNO
            
            `data` - даты наблюдений
        '''
        time = []
        date = []
        coords = []
        velocities = []
        lines = []
        megno = []
        mean_megno = []
        with open(path, 'r') as data:
            for line in data:
                line = line.strip().split(' ')
                line_clear = [x for x in line if x != '']
                lines.append(line_clear)
            
            for i in range(0, len(lines), 3):
                line1 = lines[i]
                line2 = lines[i+1]
                line3 = lines[i+2]
            
                time.append(float(line1[1])/(86400*365))
                date.append((int(line1[3]), int(line1[4]), int(line1[5])))
                
                x = float(line2[1])
                y = float(line2[2])
                z = float(line2[3])
                megno.append(float(line2[4]))

                Vx = float(line3[0])
                Vy = float(line3[1])
                Vz = float(line3[2])
                mean_megno.append(float(line3[3]))

                coords.append(np.array([x, y ,z]))
                velocities.append((Vx, Vy, Vz))
        return np.array(coords)


class Grapher:
    '''
    Класс для построения графиков
    '''
    def __init__(self, custom_rcParams=None, projection='3d'):
        if custom_rcParams:
            plt.rcParams.update(custom_rcParams)
        self.fig = plt.figure(figsize=(8,6), dpi=80)
        self.projection = projection
        self.ax = self.fig.add_subplot(1, 1, 1, projection=projection)

    def print(self, data, **kwargs):
        if self.projection == '3d':
            try:
                self.ax.plot(data[0, :], data[1, :], data[2, :], **kwargs)
            except IndexError:
                self.ax.scatter(data[0], data[1], data[2], **kwargs)

            self.ax.set_xlabel('x, км')
            self.ax.set_ylabel('y, км')
            self.ax.set_zlabel('z, км')
            self.ax.set_xlim(-40000, 40000)
            self.ax.set_ylim(-40000, 40000)
            self.ax.set_zlim(-40000, 40000)
            self.ax.set_aspect('equal', adjustable='box')
        else: 
            lmd, phi = zip(*data)
            lmd = list(map(lambda x: x.decimal, lmd))
            phi = list(map(lambda x: x.decimal, phi))

            self.ax.scatter(lmd, phi, **marker)
            self.ax.set_xlabel('$\lambda, °$')
            self.ax.set_ylabel('$\phi, °$')
        self.fig.suptitle(kwargs.get('title', ''))

    @staticmethod
    def __update(N, data, line, point):
        point._offsets3d = (data[0, N:(N+1)], data[1, N:(N+1)], data[2, N:(N+1)])
        line.set_data(data[:2, :N])
        line.set_3d_properties(data[2, :N])

    def animation(self, satellite, earth, **kwargs):
        data = np.array(list(satellite.evolution())).T
       
        self.ax.plot(*zip(*earth.coords))
        point = self.ax.scatter(data[0, 0:1], data[1, 0:1], data[2, 0:1], c='gray')
        line, = self.ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], color='k') 
        
        self.ax.set_xlim3d([-40000.0, 40000.0])
        self.ax.set_xlabel('x. км')

        self.ax.set_ylim3d([-40000.0, 40000.0])
        self.ax.set_ylabel('y, км')

        self.ax.set_zlim3d([-40000.0, 40000.0])
        self.ax.set_zlabel('z, км')
        self.ax.set_aspect('equal', adjustable='box')

        self.fig.suptitle(satellite.type + ' спутник')

        N = 100
        ani = FuncAnimation(self.fig, self.__update, N, fargs=(data, line, point), interval=1500/N, blit=False)
        if kwargs.get('save'):
            ani.save(f'animations/{kwargs.get("title", satellite.type)}-anim.gif', writer='imagemagick')
        plt.show()


    def show(self):
        plt.show()
