import numpy as np
import scipy as sp
import scipy.sparse


class StencilArray(object):

    def __init__(self):
        self.vals = np.zeros(9)

    @property
    def C(self):
        return self.vals[0]

    @C.setter
    def C(self, val):
        self.vals[0] = val

    @property
    def N(self):
        return self.vals[1]

    @N.setter
    def N(self, val):
        self.vals[1] = val

    @property
    def S(self):
        return self.vals[2]

    @S.setter
    def S(self, val):
        self.vals[2] = val

    @property
    def E(self):
        return self.vals[3]

    @E.setter
    def E(self, val):
        self.vals[3] = val

    @property
    def W(self):
        return self.vals[4]

    @W.setter
    def W(self, val):
        self.vals[4] = val

    @property
    def NE(self):
        return self.vals[5]

    @NE.setter
    def NE(self, val):
        self.vals[5] = val

    @property
    def NW(self):
        return self.vals[6]

    @NW.setter
    def NW(self, val):
        self.vals[6] = val

    @property
    def SE(self):
        return self.vals[7]

    @SE.setter
    def SE(self, val):
        self.vals[7] = val

    @property
    def SW(self):
        return self.vals[8]

    @SW.setter
    def SW(self, val):
        self.vals[8] = val


class DiffOp:

    def __init__(self, nx, ny):
        self.iv = []
        self.jv = []
        self.dv = []
        self.nx = nx
        self.ny = ny
        self.a = None

    def insert(self, i, j, stencil_array):
        nx = self.nx
        ny = self.ny
        sta = stencil_array
        tol = 1e-10

        self.iv.append(j*nx + i)
        self.jv.append(j*nx + i)
        self.dv.append(sta.C)

        if j + 1 < ny:
            if abs(sta.N) > tol:
                self.iv.append(j*nx + i)
                self.jv.append((j+1)*nx + i)
                self.dv.append(sta.N)
            if abs(sta.NE) > tol and i+1 < nx:
                self.iv.append(j*nx + i)
                self.jv.append((j+1)*nx + i + 1)
                self.dv.append(sta.NE)
            if abs(sta.NW) > tol and i > 0:
                self.iv.append(j*nx + i)
                self.jv.append((j+1)*nx + i-1)
                self.dv.append(sta.NW)
        if j > 0:
            if abs(sta.S) > tol:
                self.iv.append(j*nx + i)
                self.jv.append((j-1)*nx + i)
                self.dv.append(sta.S)
            if abs(sta.SE) > tol and i+1 < nx:
                self.iv.append(j*nx + i)
                self.jv.append((j-1)*nx + i + 1)
                self.dv.append(sta.SE)
            if abs(sta.SW) > tol and i > 0:
                self.iv.append(j*nx + i)
                self.jv.append((j-1)*nx + i - 1)
                self.dv.append(sta.SW)
        if i + 1 < nx:
            if abs(sta.E) > tol:
                self.iv.append(j*nx + i)
                self.jv.append(j*nx + i + 1)
                self.dv.append(sta.E)
        if i > 0:
            if abs(sta.W) > tol:
                self.iv.append(j*nx + i)
                self.jv.append(j*nx + i - 1)
                self.dv.append(sta.W)

    def assemble(self):
        a = sp.sparse.coo_matrix((self.dv, (self.iv, self.jv)),
                                 shape=(self.nx*self.ny, self.nx*self.ny))
        self.a = a.tocsr()


def create_op(met, eps):
    nx = met.x.shape[1]
    ny = met.x.shape[0]

    deta = 2. / (np.sqrt(2))
    ideta2 = 1. / (deta**2)
    op = DiffOp(nx-2, ny-2)
    for j in range(ny - 2):
        jj = j+1
        for i in range(nx - 2):
            ii = i+1
            sta = StencilArray()
            sta.N = met.yhJ[j+1, i] * met.g22[j+1, i]
            sta.S = met.yhJ[j, i] * met.g22[j, i]
            sta.E = met.xhJ[j, i+1] * met.g11[j, i+1]
            sta.W = met.xhJ[j, i] * met.g11[j, i]
            sta.NE = .25 * ((met.cJ[jj, ii+1] * met.g12[jj, ii+1]) +
                            (met.cJ[jj+1, ii] * met.g21[jj+1, ii]))
            sta.SE = .25 * ((-1 * met.cJ[jj, ii+1] * met.g12[jj, ii+1]) +
                            (-1 * met.cJ[jj-1, ii] * met.g21[jj-1, ii]))
            sta.NW = .25 * ((-1 * met.cJ[jj, ii-1] * met.g12[jj, ii-1]) +
                            (-1 * met.cJ[jj+1, ii] * met.g21[jj+1, ii]))
            sta.SW = .25 * ((met.cJ[jj, ii-1] * met.g12[jj, ii-1]) +
                            (met.cJ[jj-1, ii] * met.g21[jj-1, ii]))
            sta.C = ((-1. * met.xhJ[j, i+1] * met.g11[j, i+1]) +
                     (-1. * met.xhJ[j, i] * met.g11[j, i]) +
                     (-1. * met.yhJ[j+1, i] * met.g22[j+1, i]) +
                     (-1. * met.yhJ[j, i] * met.g22[j, i]))

            sta.vals = -1*sta.vals

            op.insert(i, j, sta)

    op.assemble()
    return op
