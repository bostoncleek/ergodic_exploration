import numpy as np
import matplotlib.pyplot as plt

from barrier import Barrier

class ErgodicControlKL(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples):
        self.explr_space = explr_space
        self.model = model
        self.target_dist = target_dist
        self.horizon = horizon
        self.dt = 0.1
        self.N = int(self.horizon / self.dt)
        self.num_samples = num_samples

        self.barrier = Barrier(explr_space)

        # Kinematic cart the stabilizing policy is all zeros
        self.mux = np.zeros(model.action_space_dim)
        self.u_seq = np.zeros((model.action_space_dim, self.N))
        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), 0.0)))

        self.R = np.array([[0.01, 0.0],
                           [0.0, 0.001]])

        self.q = 1.0

        self.Rinv = np.linalg.inv(self.R)

        self.pT = np.zeros(model.state_space_dim)


    def ergodic_cost_deriv(self, s, ps, x):
        sum = np.zeros(self.model.state_space_dim)
        for i in range(self.num_samples):
            sum[self.model.explr_dim]  = sum[self.model.explr_dim] + 2.0 * ps[i] * (s[:,i] - x[self.model.explr_dim]).dot(self.model.sigmaInv)
        return sum


    def controls(self, x):

        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(self.model.action_space_dim)

        # forward simulate trajectory
        xt = np.zeros((self.model.state_space_dim, self.N))
        fdx = []
        fdu = []
        dbar= []
        for i in range(self.N):
            xt[:,i] = x
            fdx.append(self.model.fdx(x, self.u_seq[:,i]))
            fdu.append(self.model.fdu(x))
            dbar.append(self.barrier.dx(x))
            x = self.model.step(x, self.u_seq[:,i], self.dt)

        # sample points in (x,y) space
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        # plt.figure(dpi=110,facecolor='w')
        # plt.plot(xt[0], xt[1])
        # plt.show()

        # xy, vals = self.target_dist.sample_grid_spec(s.T, ps)
        # plt.contourf(*xy, vals, levels=10)
        # plt.show()


        # backwards pass
        rho = self.pT
        for i in reversed(range(self.N)):
            edx = np.zeros(self.model.state_space_dim)
            edx = self.q * self.ergodic_cost_deriv(s, ps, xt[:,i])

            bdx = np.zeros(self.model.state_space_dim)
            bdx[self.model.explr_dim] = dbar[i]
            # print("edx: ", edx, " bdx: ", bdx)

            # ergodic metric is not used to update heading
            # rho = rho - self.dt * (-edx - np.dot(fdx[i].T, rho))
            rho = rho - self.dt * (-edx - bdx - np.dot(fdx[i].T, rho))

            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

        # print(self.u_seq)

        return self.u_seq[:,0].copy()
