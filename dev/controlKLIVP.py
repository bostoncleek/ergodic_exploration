import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from barrier import Barrier

class ErgodicControlKLIVP(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples, buffer_size):
        self.explr_space = explr_space
        self.model = model
        self.target_dist = target_dist
        self.horizon = horizon
        self.dt = 0.1
        self.N = int(self.horizon / self.dt)
        self.tspan = [0.0, horizon]
        self.tvec = np.linspace(0.0, horizon, self.N)
        self.num_samples = num_samples
        self.buffer_size = buffer_size
        self.barrier = Barrier(explr_space)
        self.u_seq = np.zeros((model.action_space_dim, self.N))

        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), -0.1)))

        # Penalty on controls
        # self.R = np.array([[0.1, 0.0],
        #                    [0.0, 0.01]])
        self.R = np.eye(4) * 0.01

        self.Rinv = np.linalg.inv(self.R)

        # uncertainty in (x,y) position
        self.sigma = np.eye(2) * 0.01
        self.sigmaInv = np.linalg.inv(self.sigma)

        self.pT = np.zeros(model.state_space_dim)

        self.past_states = None
        self.stored = 0


    def xdot(self, t, x):
        i = int(round((self.N-1)*(t-self.tspan[0])/(self.tspan[1]-self.tspan[0])))
        if (i == self.N or i > self.N):
            print("WARNING i is out of bounds")
            i = self.N-1
        return self.model.f(x, self.u_seq[:,i])


    def rhodot(self, t, rho, xt, s, ps, q):
        i = int(round((self.N-1)*(t-self.tspan[0])/(self.tspan[1]-self.tspan[0])))
        if (i == self.N or i > self.N):
            print("WARNING i is out of bounds")

        edx = np.zeros(self.model.state_space_dim)
        for j in range(self.num_samples):
            dq = self.traj_stat_deriv(s[:,j], xt[:,i])
            edx += (ps[j]/q[j]) * dq

        bdx = np.zeros(self.model.state_space_dim)
        bdx[self.model.explr_dim] = self.barrier.dx(xt[:,i])

        fdx = self.model.fdx(xt[:,i], self.u_seq[:,i])

        return edx - bdx - np.dot(fdx.T, rho)


    def traj_stat(self, si, xt):
        """
        Compose time average trajectory statistics
        Args:
            si: point (x,y) uniformly sampled from exploration domain
            xt: trajectory
        Returns:
            q: time average trajectory statistics (Eqn 3)
        """
        q = 0.0
        for i in range(xt.shape[1]):
            # NOTE: self.model.explr_dim = [0, 1], the exploration domain
            # is less than the dimension of the robot's state space.
            # The robot states space inludes heading but the ergodic measure does not
            # consider the heading.
            q += np.exp(-0.5 * np.dot(si - xt[self.model.explr_dim, i], self.sigmaInv).dot(si - xt[self.model.explr_dim, i]))
        return q #* 1./(xt.shape[1])


    def traj_stat_deriv(self, si, x):
        """
        Compose derivative of the time average trajectory statistics wrt the state
        Args:
            si: point (x,y) uniformly sampled from exploration domain
            x: state
        Returns:
            qx: derivative of trajectory statistics
        """
        dq = np.zeros(self.model.state_space_dim)
        q_innerds = np.exp(-0.5 * np.dot(si - x[self.model.explr_dim], self.sigmaInv).dot(si - x[self.model.explr_dim]))
        dq[self.model.explr_dim] = q_innerds * np.dot(si - x[self.model.explr_dim], self.sigmaInv)
        return dq


    def controls(self, x):
        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(self.model.action_space_dim)

        # forward simulate trajectory
        xt = solve_ivp(self.xdot, self.tspan, x, t_eval=self.tvec).y

        total_traj = None
        if self.past_states is not None:
            total_traj = np.append(self.past_states.copy(), xt, axis=1)
        else:
            total_traj = xt

        # sample points in (x,y) space
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        q = np.zeros(self.num_samples)
        for j in range(self.num_samples):
            q[j] = self.traj_stat(s[:,j], total_traj)
        q /= np.sum(q)

        # backwards pass
        rho = solve_ivp(self.rhodot, self.tspan[::-1], self.pT, t_eval=self.tvec[::-1], args=(xt, s, ps, q)).y
        rho = np.flip(rho, axis=1)

        for i in range(self.N):
            fdu = self.model.fdu(xt[:,i])

            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu.T), rho[:,i])

            if (np.abs(self.u_seq[:,i]) > 1.0).any() or (np.abs(self.u_seq[:,i]) < -1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

        if self.past_states is None:
            self.past_states = x.reshape(self.model.state_space_dim,1)
        elif self.stored < self.buffer_size :
            self.past_states = np.append(self.past_states, x.reshape(self.model.state_space_dim,1), axis=1)
            self.stored += 1
        else:
            self.past_states[:,:-1] = self.past_states[:,1:]
            self.past_states[:,-1] = x

        # print(self.u_seq)
        return self.u_seq[:,0].copy()
