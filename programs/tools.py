import matplotlib.pyplot as plt

from config import * 


class Earth:
    Mass = MASS
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
                data.append(line_filtered)
        return data

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
    def __init__(self, custom_rcParams=None):
        if custom_rcParams:
            plt.rcParams.update(custom_rcParams)
        self.fig = plt.figure(figsize=(8,6), dpi=80)
        self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')

    def print(self, data, **kwargs):
        try:
            x, y, z = zip(*data)
            self.ax.plot(x, y, z, **kwargs)
        except TypeError:
            x, y, z = data
            self.ax.scatter(x, y, z,**kwargs)
    
        self.fig.suptitle(kwargs.get('title', ''))
        self.ax.set_xlabel('x, км')
        self.ax.set_ylabel('y, км')
        self.ax.set_zlabel('z, км')
        self.ax.set_xlim(-40000, 40000)
        self.ax.set_ylim(-40000, 40000)
        self.ax.set_zlim(-40000, 40000)
        self.ax.set_aspect('equal', adjustable='box')

    def show(self):
        plt.show()
