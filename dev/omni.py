import numpy as np

class Omni(object):
    def __init__(self):
        self.action_space_dim = 4
        self.state_space_dim = 3
        # state index 0 and 1 -> (x, y) these are the states fro exploring
        self.explr_dim = np.array([0, 1])

        self.r = 0.1
        self.lx = 0.5
        self.ly = 0.5


    def fdx(self, x, u):
        s_yth = self.r/4.0 * np.sin(x[1]*x[2])
        c_yth = self.r/4.0 * np.cos(x[1]*x[2])

        df0dth = u[0]*(-s_yth + c_yth) + u[1]*(-s_yth - c_yth) + u[2]*(-s_yth - c_yth) + u[3]*(-s_yth + c_yth)

        df1dth = u[0]*(-s_yth - c_yth) + u[1]*(s_yth - c_yth) + u[2]*(s_yth - c_yth) + u[3]*(-s_yth - c_yth)

        A = np.array([[0.0, 0.0, df0dth],
                      [0.0, 0.0, df1dth],
                      [0.0, 0.0, 0.0]])

        return A


    def fdu(self, x):
        s_yth = self.r/4.0 * np.sin(x[1]*x[2])
        c_yth = self.r/4.0 * np.cos(x[1]*x[2])

        l = self.r/(4.0 * (self.lx + self.ly))

        B = np.array([[s_yth+c_yth, -s_yth+c_yth, -s_yth+c_yth, s_yth+c_yth],
                      [-s_yth+c_yth, -s_yth-c_yth, -s_yth-c_yth, -s_yth+c_yth],
                      [-l,            l,            -l,           l]])
        return B


    def f(self, x, u):
        xdot = np.zeros(3)
        s_yth = np.sin(x[1]*x[2])
        c_yth = np.cos(x[1]*x[2])

        s_yth_u0 = u[0]*s_yth
        c_yth_u0 = u[0]*c_yth

        s_yth_u1 = u[1]*s_yth
        c_yth_u1 = u[1]*c_yth

        s_yth_u2 = u[2]*s_yth
        c_yth_u2 = u[2]*c_yth

        s_yth_u3 = u[3]*s_yth
        c_yth_u3 = u[3]*c_yth

        l = 1.0/(self.lx + self.ly)

        xdot[0] = s_yth_u0 + c_yth_u0 - s_yth_u1 + c_yth_u1 - s_yth_u2 + c_yth_u2 + s_yth_u3 + c_yth_u3

        xdot[1] = -s_yth_u0 + c_yth_u0 - s_yth_u1 - c_yth_u1 - s_yth_u2 - c_yth_u2 - s_yth_u3 + c_yth_u3

        xdot[2] = -u[0]*l + u[1]*l - u[2]*l + u[3]*l

        # Rinv = np.array([[np.cos(x[2]),  np.sin(x[2]), 0.0],
        #                  [-np.sin(x[2]), np.cos(x[2]), 0.0],
        #                  [0.0,           0.0,          1.0]])
        #
        # var = 1.0/(self.lx + self.ly)
        #
        # Hp = (self.r/4.0) * np.array([[1.0,  1.0,  1.0, 1.0],
        #                               [1.0, -1.0, -1.0, 1.0],
        #                               [-var, var, -var, var]])
        #
        # return np.dot(Rinv, Hp).dot(u)

        return self.r/4.0 * xdot

    def step(self, x, u, dt):
        x = x + self.f(x, u) * dt
        return x
