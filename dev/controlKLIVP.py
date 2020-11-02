import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from barrier import Barrier

class ErgodicControlKLIVP(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples):
        self.explr_space = explr_space
        self.model = model
        self.target_dist = target_dist
        self.horizon = horizon
        self.dt = 0.01
        self.N = int(self.horizon / self.dt)
        self.tspan = [0.0, horizon]
        self.tvec = np.linspace(0.0, horizon, self.N)

        self.num_samples = num_samples

        self.barrier = Barrier(explr_space)

        self.u_seq = np.zeros((model.action_space_dim, self.N))

        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), -0.1)))

        # Penalty on controls
        self.R = np.array([[0.5, 0.0],
                           [0.0, 0.01]])

        self.q = 1.0

        self.Rinv = np.linalg.inv(self.R)

        # uncertainty in (x,y) position
        self.sigma = np.eye(2) * 0.1

        self.sigmaInv = np.linalg.inv(self.sigma)

        # self.eta = np.linalg.det(2.0*np.pi*self.sigma)**0.5

        self.pT = np.zeros(model.state_space_dim)


    # def ergodic_measure_terms(self, s, ps, x):
    #     # pdot based on eqn 20
    #     sum = np.zeros(self.model.state_space_dim)
    #     for i in range(self.num_samples):
    #         sum[self.model.explr_dim]  +=  2.0 * ps[i] * np.dot(s[:,i] - x[self.model.explr_dim], self.sigmaInv)
    #     # Not sure about the negaive sign but seems to drive it in the right direction ????
    #     return -1.*sum


    def xdot(self, t, x):
        i = int(round((self.N-1)*(t-self.tspan[0])/(self.tspan[1]-self.tspan[0])))
        if (i == self.N or i > self.N):
            print("WARNING i is out of bounds")
            i = self.N-1
        return self.model.f(x, self.u_seq[:,i])


    def rhodot(self, t, rho, xt, edx):
        i = int(round((self.N-1)*(t-self.tspan[0])/(self.tspan[1]-self.tspan[0])))
        if (i == self.N or i > self.N):
            print("WARNING i is out of bounds")

        bdx = np.zeros(self.model.state_space_dim)
        bdx[self.model.explr_dim] = self.barrier.dx(xt[:,i])

        fdx = self.model.fdx(xt[:,i], self.u_seq[:,i])

        return -edx - bdx - np.dot(fdx.T, rho)


    def ergodic_measure_terms2(self, si, xt):
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
        xt = solve_ivp(self.xdot, self.tspan, x, t_eval=self.tvec).y


        # sample points in (x,y) space
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        # xy, vals = self.target_dist.sample_grid_spec(s.T, ps)
        # plt.contourf(*xy, vals, levels=10)
        # plt.show()

        edx = np.zeros(self.model.state_space_dim)
        for i in range(self.num_samples):
            q, dq = self.ergodic_measure_terms2(s[:,i], xt)
            edx += (ps[i]/q) * dq

        edx = self.q * edx
        # print(edx)

        # backwards pass
        rho = solve_ivp(self.rhodot, self.tspan[::-1], self.pT, t_eval=self.tvec[::-1], args=(xt, edx)).y
        rho = np.flip(rho, axis=1)

        for i in range(self.N):
            fdu = self.model.fdu(xt[:,i])

            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu.T), rho[:,i])

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])


        # print(self.u_seq)
        return self.u_seq[:,0].copy()
