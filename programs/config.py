from AnglesPy import Angles

date = '02.03.2024 11:00'

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

# graph parameters
custom_rcParams = {
    'font.size': 10,
    'lines.linewidth': 0.7,
    'lines.color': 'black',
    'font.family': 'Times New Roman',
    'mathtext.fontset': 'custom',
    'mathtext.it': 'Times New Roman:italic',
}

marker = {
    'c': 'k',
    's': 50,
    'marker': '.'
}

# file paths
EARTH_PATH = 'files/Шар_м.DAT'
COEFS = 'files/garm60.txt'