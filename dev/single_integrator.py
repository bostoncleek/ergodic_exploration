import numpy as np

class SingleIntegrator(object):
    def __init__(self):
        # State linearization
        self.A = np.zeros((2,2))
        # Control linearization
        self.B = np.eye(2)

    def fdx(self):
        return self.A.copy()

    def fdu(self):
        return self.B.copy()

    def f(self, t, x, u, t0, tf, N):
        i = int(round((N-1)*(t-t0)/(tf-t0)))
        return np.dot(self.A, x) + np.dot(self.B, u[:,i])

    def step(self, x, u, dt):
        return x + (np.dot(self.A, x) + np.dot(self.B, u)) * dt
