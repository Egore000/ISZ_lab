from AnglesPy import Angles

date = '23.02.2024 18:30'

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

# accuracy
eps = 1e-13
t0 = 0

# constants
MU = 3.986004418E5
G = 6.67e-17            # H * km**2 / kg**2
MASS = 2e24             # kg
T0 = 2456514.5          # s

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