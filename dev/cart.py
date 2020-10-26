import numpy as np

class Cart(object):
    def __init__(self):
        self.action_space_dim = 2
        self.state_space_dim = 3
        # state index 0 and 1 -> (x, y) these are the states fro exploring
        self.explr_dim = np.array([0, 1])

    def fdx(self, x, u):
        A =  np.array([
            [ 0., 0., -np.sin(x[2])*u[0]],
            [ 0., 0.,  np.cos(x[2])*u[0]],
            [ 0., 0.,                0.]])
        return A

    def fdu(self, x):
        B = np.array([
            [np.cos(x[2]), 0.],
            [np.sin(x[2]), 0.],
            [0., 1.]])
        return B

    def f(self, x, u):
        return np.array([np.cos(x[2])*u[0], np.sin(x[2])*u[0], u[1]])

    def step(self, x, u, dt):
        x = x + self.f(x, u) * dt
        # x[0] += np.random.normal(0.0, 0.1**2, 1)
        # x[1] += np.random.normal(0.0, 0.1**2, 1)
        # x[2] += np.random.normal(0.0, 0.001**2, 1)
        # x[2] = x[2] % (2.0*np.pi)
        return x
