import numpy as np
import scipy as sp
import scipy.sparse as sparse
from pyamg import smoothed_aggregation_solver
from pideal import pideal_solver
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d.axes3d import Axes3D

class MyAxes3D(Axes3D):
    def __init__(self, baseObject, sides_to_draw):
        self.__class__ = type(baseObject.__class__.__name__,
                              (self.__class__, baseObject.__class__),{})
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
        
        if 'l' in self.sides_to_draw :
            # draw zaxis on the left side
            zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
                             tmp_planes[0], tmp_planes[1],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        if 'r' in self.sides_to_draw :
            # draw zaxis on the right side
            zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                             tmp_planes[1], tmp_planes[0],
                             tmp_planes[4], tmp_planes[5])
            zaxis.draw(renderer)
        
        zaxis._PLANES = tmp_planes
        
        # disable draw grid
        zaxis.axes._draw_grid = draw_grid_old


# Create mesh
n = 30;
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
X, Y = np.meshgrid(x, y)

# A = ...
ml = smoothed_aggregation_solver(A)

# Select basis function
basis = 59
P = ml.levels[0].P.todense()
Pj = P[:,basis].reshape((n,n))

plt.pcolor(xv[1:-1,1:-1], yv[1:-1,1:-1], Pj)
# Main graph
# fig = plt.figure(figsize=(12,9))
# ax = fig.add_subplot(1, 1, 1, projection='3d')
# p = ax.plot_surface(X, Y, Pj,
#                     rstride=1, cstride=1, vmin=0.0, vmax=1.0, linewidth=1,
#                     cmap=plt.cm.rainbow, alpha=0.9,  shade=True)
# plt.colorbar(p, use_gridspec=False, ax=ax)
# ax.set_zlim([0.0, 1.0])
# ax.view_init(elev=20., azim=130)
# ax._axis3don = True
# ax = fig.add_axes(MyAxes3D(ax, 'l'))
# fig.savefig(dir + solver + '_d' + str(degree) + '_j' + str(basis) + '.png',
#             bbox_inches='tight', dpi=200, transparent=True, pad_inches=0)
