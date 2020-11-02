import numpy as np

# Penalize trajectories that leave the exploration domain
class Barrier(object):
    def __init__(self, explr_space, pow=2, weight=100.0):
        self.explr_space = explr_space
        self.dl  = explr_space[1] - explr_space[0]
        self.pow = pow
        self.weight = weight
        self.eps = 0.01

    def cost(self, state):
        """
        Returns the actual cost of the barrier

        Note: only consider (x,y) of state
        """
        x = np.array([state[0], state[1]])
        cost = 0.
        cost += np.sum((x > self.explr_space[1] - self.eps) * (x - (self.explr_space[1] - self.eps))**self.pow)
        cost += np.sum((x < self.explr_space[0] + self.eps) * (x - (self.explr_space[0] + self.eps))**self.pow)
        return self.weight * cost

    def dx(self, state):
        """
        Returns the derivative of the barrier wrt to the exploration
        state

        Note: only consider (x,y) of state
        """
        x = np.array([state[0], state[1]])

        # if((x > self.explr_space[1]).any() or (x < self.explr_space[0]).any()):
        #     print("Out of bounds")

        dx = np.zeros(x.shape)
        dx += 2*(x > (self.explr_space[1] - self.eps)) * (x - (self.explr_space[1] - self.eps))
        dx += 2*(x < (self.explr_space[0] + self.eps)) * (x - (self.explr_space[0] + self.eps))
        return self.weight * dx
