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

        self.u_seq = np.zeros((model.action_space_dim, self.N))
        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), 0.0)))

        # self.R = np.array([[0.5, 0.0],
        #                    [0.0, 0.01]])

        self.R = np.array([[0.01, 0.0],
                           [0.0, 0.001]])

        self.q = 1.0

        self.Rinv = np.linalg.inv(self.R)

        # uncertainty in (x,y) position
        # self.sigma = np.eye(2)
        self.sigma = np.eye(2) * 0.1

        self.sigmaInv = np.linalg.inv(self.sigma)

        # self.eta = np.linalg.det(2.0*np.pi*self.sigma)**0.5

        self.pT = np.zeros(model.state_space_dim)


    def ergodic_cost_deriv(self, s, ps, x):
        # pdot based on eqn 20
        sum = np.zeros(self.model.state_space_dim)
        for i in range(self.num_samples):
            sum[self.model.explr_dim]  +=  2.0 * ps[i] * np.dot(s[:,i] - x[self.model.explr_dim], self.sigmaInv)
        # Not sure about the negaive sign but seems to drive it in the right direction ????
        return -1.*sum


    def ergodic_cost_deriv_long(self, si, xt):
        # pdot based on eqn 10
        # q based on Eqn 3
        q = 0.0
        dq = np.zeros(self.model.state_space_dim)
        for i in range(self.N):
            q_integrand = np.exp(-0.5 * np.dot(si - xt[self.model.explr_dim, i], self.sigmaInv).dot(si - xt[self.model.explr_dim, i]))
            dq_integrand = q_integrand * np.dot(si - xt[self.model.explr_dim, i], self.sigmaInv)

            q += q_integrand
            dq[self.model.explr_dim] += dq_integrand

        # q = (1./(self.N * self.eta)) * q
        # dq = (1./(self.N * self.eta)) * dq

        return q, -1.*dq


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

        # xy, vals = self.target_dist.sample_grid_spec(s.T, ps)
        # plt.contourf(*xy, vals, levels=10)
        # plt.show()


        edx = np.zeros(self.model.state_space_dim)
        for i in range(self.num_samples):
            q, dq = self.ergodic_cost_deriv_long(s[:,i], xt)
            edx += (ps[i]/q ) * dq

        edx = self.q * edx
        # print(edx)

        # backwards pass
        rho = self.pT
        for i in reversed(range(self.N)):
            # edx = np.zeros(self.model.state_space_dim)
            # edx = self.q * self.ergodic_cost_deriv(s, ps, xt[:,i])

            bdx = np.zeros(self.model.state_space_dim)
            bdx[self.model.explr_dim] = dbar[i]

            # ergodic metric is not used to update heading
            ## SIGN ERROR  on edx ?????
            rho = rho - self.dt * (edx -bdx -np.dot(fdx[i].T, rho))

            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

        # print(self.u_seq)

        return self.u_seq[:,0].copy()
