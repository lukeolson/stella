from cedar.cdr2 import *
from diffop import StencilArray, DiffOp

class CedarMat(DiffOp):
    """
    Attributes
    ----------
    nx, ny : int
        dimensions of the grid

    A : cedar structured matrix

    Methods
    -------
    insert(i, j, stencil_array)
        insert a stencil array at i, j

    Examples
    --------
    >>> op = DiffOp(nx-2, ny-2)
    """

    def __init__(self, nx, ny):
        """
        nx, ny : int
            grid dimensions
        """
        super(CedarMat, self).__init__(nx, ny)
        self.A = stencil_op_nine(nx, ny)
        self.A.toarray()[:,:,:] = np.zeros((nx+2,ny+2,5))


    def insert(self, i, j, stencil_array):
        """
        i, j : int
            location

        stencil_array : StencilArray
            stencil
        """

        ii = i+1
        jj = j+1
        nx = self.nx
        ny = self.ny
        sta = stencil_array
        tol = 1e-10

        # get numpy array
        op = self.A.toarray()

        op[ii,jj,nine_pt.c] = sta.C

        if j < ny:
            if abs(sta.N) > tol:
                op[ii, jj+1, nine_pt.s] = abs(sta.N)
            if abs(sta.NE) > tol and i+1 < nx:
                op[ii+1, jj+1, nine_pt.sw] = abs(sta.NE)
            if abs(sta.NW) > tol and i > 0:
                op[ii, jj+1, nine_pt.nw] = abs(sta.NW)
        if j > 1:
            if abs(sta.S) > tol:
                op[ii, jj, nine_pt.s] = abs(sta.S)
            if abs(sta.SE) > tol and i+1 < nx:
                op[ii+1, jj, nine_pt.nw] = abs(sta.SE)
            if abs(sta.SW) > tol and i > 0:
                op[ii, jj, nine_pt.sw] = abs(sta.SW)
        if i < nx:
            if abs(sta.E) > tol:
                op[ii+1, jj, nine_pt.w] = abs(sta.E)
        if i > 1:
            if abs(sta.W) > tol:
                op[ii, jj, nine_pt.w] = abs(sta.W)
