import numpy as np

class Cart(object):
    def __init__(self):
        self.action_space_dim = 2
        self.state_space_dim = 3
        # state index 0 and 1 -> (x, y) these are the states for exploring
        self.explr_dim = np.array([0, 1])

        self.wheel_radius = 0.033
        self.wheel_base = 0.08

        # self.wheel_radius = 0.2
        # self.wheel_base = 0.5

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

    # def fdx(self, x, u):
    #     A =  np.array([
    #         [ 0., 0., -(self.wheel_radius / 2.0) * (u[0] + u[1]) * np.sin(x[2])],
    #         [ 0., 0.,  (self.wheel_radius / 2.0) * (u[0] + u[1]) * np.cos(x[2])],
    #         [ 0., 0.,                0.]])
    #     return A
    #
    # def fdu(self, x):
    #     B = np.array([
    #         [np.cos(x[2]), np.cos(x[2])],
    #         [np.sin(x[2]), np.sin(x[2])],
    #         [-1.0 /self. wheel_base, 1.0 / self.wheel_base]])
    #     return (self.wheel_radius / 2.0) * B
    #
    # def f(self, x, u):
    #     xdot = np.zeros(3)
    #     xdot[0] = (u[0] + u[1]) * np.cos(x[2]);
    #     xdot[1] = (u[0] + u[1]) * np.sin(x[2]);
    #     xdot[2] = (u[1] - u[0]) / self.wheel_base;
    #     return (self.wheel_radius / 2.0) * xdot

    def step(self, x, u, dt):
        x = x + self.f(x, u) * dt
        return x



# model = Cart()
# x = np.array([1.0, 2.0, 0.707])
# u = np.array([0.5, 0.01])
#
# print(model.f(x,u))
# print(model.fdx(x,u))
# print(model.fdu(x))
