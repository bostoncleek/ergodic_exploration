import numpy as np

class SingleIntegrator(object):
    def __init__(self):
        self.action_space_dim = 2
        self.state_space_dim = 2
        # state index 0 and 1 -> (x, y) these are the states fro exploring
        self.explr_dim = np.array([0, 1])

        # uncertainty in (x,y) position
        # self.sigma = np.eye(2)
        self.sigma = np.eye(2) * 0.01
        self.sigmaInv = np.linalg.inv(self.sigma)

        # State linearization
        self.A = np.zeros((2,2))
        # Control linearization
        self.B = np.eye(2)


    def fdx(self, x, u):
        return self.A.copy()

    def fdu(self, x):
        return self.B.copy()

    def f(self, x, u):
        return np.dot(self.A, x) + np.dot(self.B, u)

    def step(self, x, u, dt):
        return x + (np.dot(self.A, x) + np.dot(self.B, u)) * dt
