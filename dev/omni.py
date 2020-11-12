import numpy as np

class Omni(object):
    def __init__(self):
        # self.action_space_dim = 4
        # self.state_space_dim = 3

        self.action_space_dim = 3
        self.state_space_dim = 3

        # state index 0 and 1 -> (x, y) these are the states fro exploring
        self.explr_dim = np.array([0, 1])

        self.r = 0.2
        self.lx = 0.5
        self.ly = 0.5

    def f(self, x, u):
        xdot = np.array([u[0]*np.cos(x[2]) - u[1]*np.sin(x[2]),
                         u[0]*np.sin(x[2]) + u[1]*np.cos(x[2]),
                         u[2]])
        return xdot

    def fdx(self, x, u):
        A =  np.array([[0.0, 0.0, -u[0]*np.sin(x[2]) - u[1]*np.cos(x[2])],
                       [0.0, 0.0, u[0]*np.cos(x[2]) - u[1]*np.sin(x[2])],
                       [0.0, 0.0, 0.0]])
        return A

    def fdu(self, x):
        B = np.array([[np.cos(x[2]), -np.sin(x[2]), 0.0],
                      [np.sin(x[2]), np.cos(x[2]), 0.0],
                      [0.0,          0.0,          1.0]])
        return B

    # def fdx(self, x, u):
    #     s = self.r/4.0 * np.sin(x[2])
    #     c = self.r/4.0 * np.cos(x[2])
    #
    #     df0dth = u[0]*(-s + c) + u[1]*(-s - c) + u[2]*(-s + c) + u[3]*(-s - c)
    #
    #     df1dth = u[0]*(s + c) + u[1]*(-s + c) + u[2]*(s + c) + u[3]*(-s + c)
    #
    #     A = np.array([[0.0, 0.0, df0dth],
    #                   [0.0, 0.0, df1dth],
    #                   [0.0, 0.0, 0.0]])
    #
    #     return A
    #
    #
    # def fdu(self, x):
    #     s = self.r/4.0 * np.sin(x[2])
    #     c = self.r/4.0 * np.cos(x[2])
    #
    #     l = self.r/(4.0 * (self.lx + self.ly))
    #
    #     B = np.array([[s+c, -s+c, s+c, -s+c],
    #                   [s-c, s+c, s-c, s+c],
    #                   [-l,   l,    l,   -l]])
    #     return B
    #
    #
    # def f(self, x, u):
    #     xdot = np.zeros(3)
    #     s = self.r/4.0 * np.sin(x[2])
    #     c = self.r/4.0 * np.cos(x[2])
    #
    #     l = self.r/(4.0 * (self.lx + self.ly))
    #
    #     xdot[0] = u[0]*(s + c) + u[1]*(-s + c) + u[2]*(s + c) + u[3]*(-s + c)
    #
    #     xdot[1] = u[0]*(s - c) + u[1]*(s + c) + u[2]*(s - c) + u[3]*(s + c)
    #
    #     xdot[2] = -u[0]*l + u[1]*l + u[2]*l - u[3]*l
    #
    #     # # This model is based on kevin lynches book
    #     # R = np.array([[np.cos(x[2]),  -np.sin(x[2]), 0.0],
    #     #               [np.sin(x[2]), np.cos(x[2]), 0.0],
    #     #               [0.0,           0.0,          1.0]])
    #     #
    #     # var = 1.0/(self.lx + self.ly)
    #     #
    #     # Hp = (self.r/4.0) * np.array([[1.0,  1.0,  1.0, 1.0],
    #     #                               [-1.0, 1.0, -1.0, 1.0],
    #     #                               [-var, var, var, -var]])
    #     #
    #     # return np.dot(Rinv, Hp).dot(u)
    #
    #     return xdot

    def step(self, x, u, dt):
        x = x + self.f(x, u) * dt
        return x

# model = Omni()
# x = np.array([0.04, 0.0, 0.0])
# u = np.array([-1.0, -1.0, -1.0, -1.0])
#
# print(model.f(x,u))
# print(model.fdx(x,u))
# print(model.fdu(x))
