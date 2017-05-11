from test import *


def run(n):
    """
    run a problem of size n
    """
    nx = n
    ny = n

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x, y)

    # xtmp = map_second_x(xv, yv)
    # yv = map_second_y(xv, yv)
    # xv = xtmp

    met = Metric(xv, yv)
    eps = np.ones((nx-1, ny-1))

    deta = 2. / (np.sqrt(2))
    ideta2 = 1. / (deta**2)
    op = DiffOp(nx-2, ny-2)
    for j in range(ny - 2):
        for i in range(nx - 2):
            sta = StencilArray()
            sta.N = ((ideta2 * eps[j+1, i+1] * met.J[j+1, i+1] * met.g12[j+1, i+1]) +
                     (ideta2 * eps[j+1, i+1] * met.J[j+1, i] * met.g21[j+1, i]))
            sta.S = ((ideta2 * eps[j, i] * met.J[j, i] * met.g12[j, i]) +
                     (ideta2 * eps[j, i+1] * met.J[j, i+1] * met.g21[j, i+1]))
            sta.E = (-1*(ideta2 * eps[j+1, i+1] * met.J[j+1, i+1] * met.g12[j+1, i+1]) +
                     -1*(ideta2 * eps[j, i+1] * met.J[j, i+1] * met.g21[j, i+1]))
            sta.W = (-1*(ideta2 * eps[j, i] * met.J[j, i] * met.g12[j, i]) +
                     -1*(ideta2 * eps[j+1, i] * met.J[j+1, i] * met.g21[j+1, i]))
            sta.NE = ideta2 * eps[j+1, i+1] * met.J[j+1, i+1] * met.g11[j+1, i+1]
            sta.NW = ideta2 * eps[j+1, i] * met.J[j+1, i] * met.g22[j+1, i]
            sta.SE = ideta2 * eps[j, i+1] * met.J[j, i+1] * met.g22[j, i+1]
            sta.SW = ideta2 * eps[j, i] * met.J[j, i] * met.g11[j, i]
            sta.C = ideta2 * ((-1 * eps[j+1, i+1] * met.J[j+1, i+1] * met.g11[j+1, i+1]) +
                              (-1 * eps[j, i] * met.J[j, i] * met.g11[j, i]) +
                              (-1 * eps[j+1, i] * met.J[j+1, i] * met.g22[j+1, i]) +
                              (-1 * eps[j, i+1] * met.J[j, i+1] * met.g22[j, i+1]))

            sta.vals = -1*sta.vals

            op.insert(i, j, sta)

    op.assemble()

    f = met.cJ * rhs(xv[1:-1, 1:-1], yv[1:-1, 1:-1])
    b = f.ravel()
    ml = pyamg.ruge_stuben_solver(op.A)
    M = ml.aspreconditioner(cycle='V')
    x, info = cg(op.A, b, tol=1e-8, maxiter=100, M=M)
    x = x.reshape(f.shape)
    xex = sol(xv[1:-1, 1:-1], yv[1:-1, 1:-1])

    diff = xex - x
    return np.max(diff)


def order():
    ns = [50, 100, 150, 200, 250, 300, 350, 400]
    norms = []
    for n in ns:
        nrm = run(n)
        norms.append(nrm)
    hs = [1./n for n in ns]
    np.savetxt('data/h.txt', hs)
    np.savetxt('data/norm.txt', norms)


def plot():
    hs = np.loadtxt('data/h.txt')
    norms = np.loadtxt('data/norm.txt')
    plt.loglog(hs, norms, label=r'Computed')
    plt.loglog(hs, hs**2, label=r'$O(h^2)$')
    plt.xlabel(r'$h$')
    plt.ylabel(r'$ || e ||_{\infty} $')
    plt.legend(loc=0)
    plt.show()


if __name__ == '__main__':
    order()
    plot()
