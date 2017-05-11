import numpy as np
from scipy import interpolate
from scipy.interpolate import interp1d


class Metric:

    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.compute()

    def boundary(self):
        x = self.x
        y = self.y

        x_s = np.zeros(x.shape)
        y_s = np.zeros(x.shape)
        x_t = np.zeros(x.shape)
        y_t = np.zeros(x.shape)

        sp = [3*.5, .5, -4*.5]
        sm = [-3*.5, 4*.5, -.5]

        x_s[:, 0] = sm[0] * x[:, 0] + sm[1] * x[:, 1] + sm[2] * x[:, 2]
        y_s[:, 0] = sm[0] * y[:, 0] + sm[1] * y[:, 1] + sm[2] * y[:, 2]
        x_s[:, -1] = sp[0] * x[:, -1] + sp[-1] * x[:, -2] + sp[-2] * x[:, -3]
        y_s[:, -1] = sp[0] * y[:, -1] + sp[-1] * y[:, -2] + sp[-2] * x[:, -3]

        x_t[0, :] = sm[0] * x[0, :] + sm[1] * x[1, :] + sm[2] * x[2, :]
        y_t[0, :] = sm[0] * y[0, :] + sm[1] * y[1, :] + sm[2] * y[2, :]
        x_t[-1, :] = sp[0] * x[-1, :] + sp[-1] * x[-2, :] + sp[-2] * x[-3, :]
        y_t[-1, :] = sp[0] * y[-1, :] + sp[-1] * y[-2, :] + sp[-2] * y[-3, :]

        return (x_s, y_s, x_t, y_t)

    def interp_diff(self, axis):
        # TODO: no reason this cant be a vector function

        x = self.x
        y = self.y

        if axis == 0:
            XHx_t = np.zeros((x.shape[0] - 2, x.shape[1] - 1))
            XHy_t = np.zeros((x.shape[0] - 2, x.shape[1] - 1))
            for j in range(1, x.shape[0] - 1):
                for i in range(1, x.shape[1]-1, 2):
                    jj = j-1
                    if i < x.shape[1] - 1:
                        f0 = interp1d(np.arange(3), x[j-1, i-1:i+2])
                        f1 = interp1d(np.arange(3), x[j+1, i-1:i+2])
                        XHx_t[jj, i-1] = (f1(.5) - f0(.5)) * .5
                        XHx_t[jj, i] = (f1(1.5) - f0(1.5)) * .5
                        f0 = interp1d(np.arange(3), y[j-1, i-1:i+2])
                        f1 = interp1d(np.arange(3), y[j+1, i-1:i+2])
                        XHy_t[jj, i-1] = (f1(.5) - f0(.5)) * .5
                        XHy_t[jj, i] = (f1(1.5) - f0(1.5)) * .5
                    elif i == x.shape[1] - 1:
                        f0 = interp1d(np.arange(3), x[j-1, i-2:i+1])
                        f1 = interp1d(np.arange(3), x[j+1, i-2:i+1])
                        XHx_t[jj, i-1] = (f1(.5) - f0(.5)) * .5
                        f0 = interp1d(np.arange(3), y[j-1, i-2:i+1])
                        f1 = interp1d(np.arange(3), y[j+1, i-2:i+1])
                        XHy_t[jj, i-1] = (f1(.5) - f0(.5)) * .5

            return XHx_t, XHy_t
        else:
            YHx_s = np.zeros((x.shape[0] - 1, x.shape[1] - 2))
            YHy_s = np.zeros((x.shape[0] - 1, x.shape[1] - 2))
            for j in range(1, x.shape[0] - 1, 2):
                for i in range(1, x.shape[1]-1):
                    ii = i-1
                    if j < x.shape[0] - 1:
                        f0 = interp1d(np.arange(3), x[j-1:j+2, i-1])
                        f1 = interp1d(np.arange(3), x[j-1:j+2, i+1])
                        YHx_s[j-1, ii] = (f1(.5) - f0(.5)) * .5
                        YHx_s[j, ii] = (f1(1.5) - f0(1.5)) * .5
                        f0 = interp1d(np.arange(3), y[j-1:j+2, i-1])
                        f1 = interp1d(np.arange(3), y[j-1:j+2, i+1])
                        YHy_s[j-1, ii] = (f1(.5) - f0(.5)) * .5
                        YHy_s[j, ii] = (f1(1.5) - f0(1.5)) * .5
                    elif j == x.shape[1] - 1:
                        f0 = interp1d(np.arange(3), x[j-2:j+1, i-1])
                        f1 = interp1d(np.arange(3), x[j-2:j+1, i+1])
                        YHx_s[j-1, ii] = (f1(.5) - f0(.5)) * .5
                        f0 = interp1d(np.arange(3), y[j-2:j+1, i-1])
                        f1 = interp1d(np.arange(3), y[j-2:j+1, i+1])
                        YHy_s[j-1, ii] = (f1(.5) - f0(.5)) * .5

            return YHx_s, YHy_s

    def compute(self):
        x = self.x
        y = self.y

        XHx_s = x[1:-1, 1:] - x[1:-1, :-1]
        XHy_s = y[1:-1, 1:] - y[1:-1, :-1]
        YHx_t = x[1:, 1:-1] - x[:-1, 1:-1]
        YHy_t = y[1:, 1:-1] - y[:-1, 1:-1]

        XHx_t, XHy_t = self.interp_diff(axis=0)
        YHx_s, YHy_s = self.interp_diff(axis=1)

        (Cx_s, Cy_s, Cx_t, Cy_t) = self.boundary()

        Cx_t[1:-1, :] = (x[2:, :] - x[:-2, :]) * .5
        Cy_t[1:-1, :] = (y[2:, :] - y[:-2, :]) * .5
        Cx_s[:, 1:-1] = (x[:, 2:] - x[:, :-2]) * .5
        Cy_s[:, 1:-1] = (y[:, 2:] - y[:, :-2]) * .5

        g_11 = YHx_s * YHx_s + YHy_s * YHy_s
        XHg_11 = XHx_s * XHx_s + XHy_s * XHy_s
        Cg_11 = Cx_s * Cx_s + Cy_s * Cy_s
        g_22 = XHx_t * XHx_t + XHy_t * XHy_t
        YHg_22 = YHx_t * YHx_t + YHy_t * YHy_t
        Cg_22 = Cx_t * Cx_t + Cy_t * Cy_t
        g_12 = Cx_s * Cx_t + Cy_s * Cy_t
        XHg_12 = XHx_s * XHx_t + XHy_s * XHy_t
        YHg_12 = YHx_s * YHx_t + YHy_s * YHy_t

        cJ = Cx_s * Cy_t - Cx_t * Cy_s
        xhJ = XHx_s * XHy_t - XHx_t * XHy_s
        yhJ = YHx_s * YHy_t - YHx_t * YHy_s

        XHg = XHg_11 * g_22 - XHg_12 * XHg_12
        YHg = g_11 * YHg_22 - YHg_12 * YHg_12
        Cg = Cg_11 * Cg_22 - g_12 * g_12

        g11 = 1./XHg * g_22
        g22 = 1./YHg * g_11
        g12 = -1./(cJ*cJ) * g_12

        self.g11 = g11
        self.g22 = g22
        self.g12 = g12
        self.g21 = g12

        self.cJ = cJ
        self.xhJ = xhJ
        self.yhJ = yhJ
