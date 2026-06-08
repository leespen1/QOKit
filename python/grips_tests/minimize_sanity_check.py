import scipy
import numpy as np


def g(y):
    return 2 * y[1] + 3 * y[2]


def make_f(n):
    def f(x):
        print("x = ", x)
        y = x[:n]
        z = y[n:]
        return np.linalg.norm(y) + np.linalg.norm(z)

    return f


x0 = np.hstack([np.array([1, 2]), np.array([1, 2])])
result = scipy.optimize.minimize(make_f(2), x0, method="COBYLA", args=())
