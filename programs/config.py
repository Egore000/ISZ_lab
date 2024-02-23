from AnglesPy import Angles

date = '22.02.2024 23:00'

# elements
a = 42169.14              # km
e = 0
i = Angles(10, 0, 0)    # grad
T = 24 * 86164/86400                 # h
Omega = Angles(0)      # grad
w = Angles(0)          # grad
M0 = Angles(0)         # grad

# # elements
# a = 42165              # km
# e = 0.001
# i = Angles(1, 0, 0)    # grad
# T = 24                 # h
# Omega = Angles(0)      # grad
# w = Angles(0)          # grad
# M0 = Angles(0)         # grad


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
    'font.size': 8,
    'lines.linewidth': 0.7,
    'lines.color': 'black',
    'font.family': 'Times New Roman'
}

# file paths
EARTH_PATH = 'files/Шар_м.DAT'
COEFS = 'files/garm60.txt'