import matplotlib as mpl
mpl.use('Qt4Agg')

import numpy as np
from cedar.cdr2 import *
from cedar_util import CedarMat
from diffop import create_op
from metric import Metric
import matplotlib.pyplot as plt

def rhs(x, y):
    return 8*(np.pi*np.pi)*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)

def sol(x, y):
    return np.sin(2*np.pi*x)*np.sin(2*np.pi*y)

def map_distort_x(x, y):
    eps = .05
    return x + eps*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def map_distort_y(x, y):
    eps = .05
    return y + eps*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def map_second_x(x, y):
    return 10*x


def map_second_y(x, y):
    return y

nx = 51
ny = 51

x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

# xtmp = map_distort_x(xv, yv)
# yv = map_distort_y(xv, yv)
# xv = xtmp

met = Metric(xv, yv)
eps = np.ones((ny-1, nx-1))
op = create_op(met, CedarMat)

with open('op.txt','w') as fh:
    fh.write(op.A.__repr__())
    fh.close()

b = grid_func(nx, ny)
x = grid_func(nx, ny)

# b.toarray()[:,:] = sol(xv, yv)
# b.toarray()[1:-1,1:-1] = met.cJ[1:-1, 1:-1] * rhs(xv[1:-1, 1:-1], yv[1:-1, 1:-1])


# solve(op.A, x, b)

# plt.pcolormesh(x.toarray()[1:-1,1:-1])
# plt.show()
