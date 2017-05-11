from test import *
from diffop import *


def run(n):
    """
    run a problem of size n
    """
    nx = n
    ny = n

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    xv, yv = np.meshgrid(x, y)

    met = Metric(xv, yv)

    op = create_op(met)

    f = met.cJ[1:-1, 1:-1] * rhs(xv[1:-1, 1:-1], yv[1:-1, 1:-1])
    b = f.ravel()

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

    return np.max(diff)


def order():
    ns = [51, 101, 201]
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
    plt.loglog(hs, norms, '-o', label=r'Computed')
    plt.loglog(hs, hs**2, label=r'$O(h^2)$')
    plt.xlabel(r'$h$')
    plt.ylabel(r'$ || e ||_{\infty} $')
    plt.legend(loc=0)
    plt.show()


if __name__ == '__main__':
    order()
    plot()
