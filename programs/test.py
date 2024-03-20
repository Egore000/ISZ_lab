import numpy as np


m = 536
N = 2067

I = np.random.random(m) # тут твои значения
sum_I = sum(I)

w = np.random.random(N)
sigma = np.random.random(m)

delta_V = 1
tau = 1
P = 1
c = 1

def A(w, I, sigma) -> complex:
    '''
    w - число
    I - np.array
    sigma - np.array
    '''
    sum = 0
    for k in range(len(I)):
        sum += I[k] / complex(tau, w - sigma[k] - delta_V*P)
    return sum / sum_I


def R(w, I, sigma) -> float:
    A_w = A(w, I, sigma)
    return 1/np.pi * ((A_w / (1 - c * tau * A_w)) * sum_I).real

R = [R(w_k, I, sigma) for w_k in w]