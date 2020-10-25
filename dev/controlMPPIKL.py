import numpy as np
from scipy.signal import savgol_filter

from barrier import Barrier

class ErgodicControlMPPIKL(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples):
        self.explr_space = explr_space
        self.model = model
        self.target_dist = target_dist
        self.horizon = horizon
        self.dt = 0.01
        self.N = int(self.horizon / self.dt)
        self.num_samples = num_samples

        # number of roll outs
        self.K = 5

        self.barrier = Barrier(explr_space)

        self.u_seq = np.zeros((model.action_space_dim, self.N))

        self.R = np.array([[0.5, 0.0],
                           [0.0, 0.01]])

        # Same dim as controls
        self.cov = np.array([0.2, 0.8])

        self.lam = 0.01
        # self.q = 1.0

    def loss(self, u, x, s, ps):
        ergodic_cost = 0.0
        for i in range(self.num_samples):
            ergodic_cost += ps[i] * np.dot(s[:,i] - x[self.model.explr_dim], self.model.sigmaInv).dot(s[:,i] - x[self.model.explr_dim])
        return ergodic_cost + u.T.dot(self.R).dot(u) + self.barrier.cost(x)


    def controls(self, x_curr):

        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(self.model.action_space_dim)

        # store perturbuations
        du1 = np.zeros((self.N,self.K))
        du2 = np.zeros((self.N,self.K))
        # record cost at each discretization for each rollout
        J = np.zeros((self.N,self.K))

        # sample points in (x,y) space
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        for k in range(self.K):
            x = x_curr

            # perturbuations
            deltau1 = np.random.normal(0.0, self.cov[0], self.N)
            deltau2 = np.random.normal(0.0, self.cov[1], self.N)
            du = np.vstack((deltau1, deltau2))

            du1[:,k] = du[0]
            du2[:,k] = du[1]

            # simulate dynamics
            for i in range(self.N):
                u_pert = self.u_seq[:,i] + du[:,i]
                x = self.model.step(x, u_pert, self.dt)

                J[i,k] = self.loss(u_pert, x, s, ps)

        # sum across each rollout backwards in time
        # consider how current actions effect future states
        J = np.cumsum(J[::-1], 0)[::-1,:]

        # update the controls
        for i in range(self.N):
            # min across all rollouts for this point in time
            beta = np.amin(J[i])
            J[i] -= beta

            # compose weights
            w = np.exp(-J[i]/self.lam) + 1e-8
            w /= np.sum(w)

            self.u_seq[0,i] += np.sum(np.dot(w, du1[i,:]))
            self.u_seq[1,i] += np.sum(np.dot(w, du2[i,:]))

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

        return self.u_seq[:,0].copy()
























#
