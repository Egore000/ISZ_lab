from AnglesPy import Angles


# Дата наблюдения
date = '14.03.2024 11:58:04'

# Параметры спутников
parameters = {
    'Тестовый': {
        'a': 42169.14,
        'e': 0,
        'i': Angles(10, 0, 0),
        'T': 86164,
        'Omega': Angles(0),
        'w': Angles(0),
        'M0': Angles(0),
    },
    'Геостационарный': {
        'a': 42165,
        'e': 0.001,
        'i': Angles(1, 0, 0),
        'T': 86400,
        'Omega': Angles(0),
        'w': Angles(0),
        'M0': Angles(0),
    }
}

# Параметры графиков
custom_rcParams = {
    'font.size': 14,
    'lines.linewidth': 0.7,
    'lines.color': 'black',
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'custom',
    'mathtext.it': 'Times New Roman:italic',
}

# Параметры точек на графиках
marker = {
    'c': 'k',
    's': 50,
    'marker': '.'
}

# Пути к файлам
EARTH_PATH = 'files/Шар_м.DAT'
COEFS = 'files/garm60.txt'
