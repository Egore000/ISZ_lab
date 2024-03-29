import os
from datetime import datetime, timezone, timedelta

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import requests
import json

from config import *


class Filer:
    '''Класс для чтения и записи в файлы'''
    @staticmethod
    def read(file: str) -> list:
        '''Чтение данных из файла `file`'''
        data = []
        with open(file, 'r') as f:
            for line in f:
                line_filtered = [float(x) for x in line.split()]
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
    def write(path: str, data: str):
        with open(path, 'w', encoding='UTF-8') as file:
            file.write(data)


class Grapher:
    '''
    Класс для построения графиков
    '''

    def __init__(self, custom_rcParams=None, projection='3d'):
        if custom_rcParams:
            plt.rcParams.update(custom_rcParams)
        self.fig = plt.figure(figsize=(8, 6), dpi=80)
        self.projection = projection
        self.ax = self.fig.add_subplot(1, 1, 1, projection=projection)

    def print(self, data, **kwargs):
        if kwargs.get('title'):
            title = kwargs.pop('title')
        else:
            title = ''

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
            try:
                lmd, phi = zip(*data)
                lmd = list(map(lambda x: x.decimal, lmd))
                phi = list(map(lambda x: x.decimal, phi))
            except TypeError:
                lmd, phi = data
                lmd = lmd.decimal
                phi = phi.decimal

            # self.ax.set_xlim(-180, 180)
            self.ax.scatter(lmd, phi, **kwargs)
            self.ax.set_xlabel('$\lambda, °$')
            self.ax.set_ylabel('$\phi, °$')
        self.fig.suptitle(title)

    @staticmethod
    def __update(N, data, line, point):
        point._offsets3d = (
            data[0, N:(N+1)], data[1, N:(N+1)], data[2, N:(N+1)])
        line.set_data(data[:2, :N])
        line.set_3d_properties(data[2, :N])

    def animation(self, satellite, earth, **kwargs):
        data = np.array(list(satellite.evolution())).T

        self.ax.plot(*zip(*earth.coords))
        point = self.ax.scatter(
            data[0, 0:1], data[1, 0:1], data[2, 0:1], c='gray')
        line, = self.ax.plot(
            data[0, 0:1], data[1, 0:1], data[2, 0:1], color='k')

        self.ax.set_xlim3d([-40000.0, 40000.0])
        self.ax.set_xlabel('x. км')

        self.ax.set_ylim3d([-40000.0, 40000.0])
        self.ax.set_ylabel('y, км')

        self.ax.set_zlim3d([-40000.0, 40000.0])
        self.ax.set_zlabel('z, км')
        self.ax.set_aspect('equal', adjustable='box')

        self.fig.suptitle(satellite.type + ' спутник')

        N = 100
        ani = FuncAnimation(self.fig, self.__update, N, fargs=(
            data, line, point), interval=1500/N, blit=False)
        if kwargs.get('save'):
            ani.save(
                f'animations/{kwargs.get("title", satellite.type)}-anim.gif', writer='imagemagick')
        plt.show()

    def show(self):
        plt.show()


class Parser:
    '''
    Класс для сбора данных о группировке ГЛОНАСС с сайта https://glonass-iac.ru/
    '''
    url = 'https://glonass-iac.ru/glonass/ephemeris/ephemeris_json.php'
    current_pos_url = 'https://glonass-iac.ru/glonass/currentPosition/getsatpos_iac.php'
    headers = {
        'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.77 Safari/537.36"
    }

    def get_response(self, url):
        '''
        Отправка запроса на сервер
        '''
        try:
            response = requests.get(url=url,
                                    headers=self.headers)
            return response
        except requests.exceptions.HTTPError as error:
            print(f'Ошибка: {error}')

    def get_response_from_file(self, path: str) -> dict:
        '''
        Получение данных из сохранённого `path.json` файла
        '''
        with open(path, 'r', encoding='UTF-8') as f:
            response = json.load(f)
        return response

    def write_to_file(self, data: dict, path: str):
        '''
        Запись `data` в `path.json` файл 
        '''
        with open(path, 'w', encoding='UTF-8') as f:
            json.dump(data, f)

    def __get_data(self, data):
        '''
        Приведение данных `data` к нужным типам
        '''
        result = {}
        for dict_ in data:
            dict_.pop('name')
            dict_.pop('color')
            dict_['Tomega'] = float(dict_['Tomega'])
            dict_['Tapp'] = float(dict_['Tapp'])
            dict_['e'] = float(dict_['e'])
            dict_['i'] = float(dict_['i'])
            dict_['Lomega'] = float(dict_['Lomega'])
            dict_['W'] = float(dict_['W'])
            dict_['deltaT2'] = float(dict_['deltaT2'])
            dict_['nl'] = int(dict_['nl'])
            dict_['deltaT'] = float(dict_['deltaT'])

            date = dict_['datetime']
            day, month, year = date.split('.')
            year = '20' + year
            dict_['datetime'] = '.'.join([day, month, year])

            result[dict_.pop('ns')] = dict_

        return result

    def parse_data(self):
        '''
        Основная функция для парсинга данных с сайта
        '''
        path = f'files/{datetime.now().date()}.json'
        if not os.path.exists(path):
            response = self.get_response(self.url)
            data = response.json()
            result = self.__get_data(data)
            self.write_to_file(result, path)
        else:
            result = self.get_response_from_file(path)

        for sat, data in result.items():
            result[sat]['i'] = Angles(decimal=result[sat]['i'])
            result[sat]['W'] = Angles(decimal=result[sat]['W'])
            result[sat]['Lomega'] = Angles(decimal=result[sat]['Lomega'])
        return result

    def get_current_position(self):
        response = self.get_response(self.current_pos_url)
        data = response.json()
        print(f'Получены данные на момент времени: {datetime.now()}')

        result = {}
        for item in data:
            if item['lon'] > 180:
                item['lon'] = item['lon'] - 360
            result[item['point']] = (item['lon'], item['lat'])
        return result


if __name__ == '__main__':
    ...
