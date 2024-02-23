import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from config import * 


class Earth:
    Mass = 2e24
    Radius = 6378
    def __init__(self, file: str):
        self.coords = Filer.read(file)


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
    def update(N, data, line, point):
        point._offsets3d = (data[0, N:(N+1)], data[1, N:(N+1)], data[2, N:(N+1)])
        line.set_data(data[:2, :N])
        line.set_3d_properties(data[2, :N])

    def animation(self, satellite, earth, **kwargs):
        data = np.array(list(satellite.evolution())).T

        self.ax.plot(*zip(*earth.coords))
        point = self.ax.scatter(data[0, 0:1], data[1, 0:1], data[2, 0:1], c='red')
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
        ani = FuncAnimation(self.fig, self.update, N, fargs=(data, line, point), interval=1500/N, blit=False)
        if kwargs.get('save'):
            ani.save(f'animations/{satellite.type}-anim.gif', writer='imagemagick')
        plt.show()


    def show(self):
        plt.show()
