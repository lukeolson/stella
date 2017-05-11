import numpy as np
import scipy as sp
import scipy.sparse
import matplotlib.pyplot as plt
from scipy.sparse.linalg import cg
import pyamg
from diffop import create_op, DiffOp, StencilArray
from metric import Metric
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d.axes3d import Axes3D


class MyAxes3D(Axes3D):
    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__), {})
        self.__dict__ = baseObject.__dict__
        self.sides_to_draw = list(sides_to_draw)
        self.mouse_init()

    def set_some_features_visibility(self, visible):
        for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
            t.set_visible(visible)
        self.w_zaxis.line.set_visible(visible)
        self.w_zaxis.pane.set_visible(visible)
        self.w_zaxis.label.set_visible(visible)

    def draw(self, renderer):
        # set visibility of some features False
        self.set_some_features_visibility(False)
        # draw the axes
        super(MyAxes3D, self).draw(renderer)
        # set visibility of some features True.
        # This could be adapted to set your features to desired visibility,
        # e.g. storing the previous values and restoring the values
        self.set_some_features_visibility(True)

        zaxis = self.zaxis
        draw_grid_old = zaxis.axes._draw_grid
        # disable draw grid
        zaxis.axes._draw_grid = False

        tmp_planes = zaxis._PLANES

        if 'l' in self.sides_to_draw:
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw:
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                             tmp_planes[1], tmp_planes[0],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)

        zaxis._PLANES = tmp_planes

        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old


def map_curl_x(x, y):
    return -1*y + 1


def map_curl_y(x, y):
    return x


def map_second_x(x, y):
    return 10*x


def map_second_y(x, y):
    return y


def map_distort_x(x, y):
    eps = .01
    return x + eps*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def map_distort_y(x, y):
    eps = .01
    return y + eps*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def rhs(x, y):
    return 8*(np.pi*np.pi)*np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def sol(x, y):
    return np.sin(2*np.pi*x)*np.sin(2*np.pi*y)


def jump_rhs(x, y):
    ret = np.zeros(x.shape)
    ret[x <= 0] = -16
    ret[x > 0] = 16
    return ret


def jump_sol(x, y):
    ret = np.zeros(x.shape)
    ret[x <= 0] = 2*(x[x <= 0]**2) + 2*x[x <= 0]
    ret[x > 0] = -4*(x[x > 0]**2) + 4*x[x > 0]
    return ret


nx = 101
ny = 101

x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)

xtmp = map_distort_x(xv, yv)
yv = map_distort_y(xv, yv)
xv = xtmp

met = Metric(xv, yv)
eps = np.ones((ny-1, nx-1))

op = create_op(met)

f = met.cJ[1:-1, 1:-1] * rhs(xv[1:-1, 1:-1], yv[1:-1, 1:-1])
b = f.ravel()
# ml = pyamg.ruge_stuben_solver(op.A,
#                               max_levels=10,
#                               max_coarse=5,
#                               keep=True)

ml = pyamg.smoothed_aggregation_solver(op.A,
                                       max_levels=10,
                                       max_coarse=5,
                                       strength='classical',
                                       keep=True)

M = ml.aspreconditioner(cycle='V')
res = []
x, info = pyamg.krylov.cg(op.A, b, maxiter=200, M=M, residuals=res)
# x = ml.solve(b, tol=1e-8, maxiter=200, residuals=res)
x = x.reshape(f.shape)
xex = sol(xv[1:-1, 1:-1], yv[1:-1, 1:-1])

diff = xex - x
print('Max Norm:', np.max(np.abs(diff)))
print('Residuals to converge:', len(res))
# plt.pcolormesh(xv)
# plt.show()
# plt.pcolormesh(yv)

# plt.show()
# plt.subplot(221)
# plt.pcolormesh(x)
# plt.title('Computed Solution')
# plt.colorbar()
# plt.subplot(222)
# plt.pcolormesh(xex)
# plt.title('Exact Solution')
# plt.colorbar()
# plt.subplot(223)
# plt.pcolormesh(diff)
# plt.title('Error')
# plt.colorbar()
# plt.show()
plt.spy(op.A - op.A.T)
plt.show()
# np.savetxt('nested.txt', res)
