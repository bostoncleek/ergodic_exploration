import numpy as np
import matplotlib.pyplot as plt

from fourier import Basis
from barrier import Barrier

class ErgodicControlSE2(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_basis):
        self.explr_space = explr_space
        self.model = model
        self.target_dist = target_dist
        self.horizon = horizon
        self.dt = 0.1
        self.N = int(self.horizon / self.dt)

        self.basis = Basis(explr_space, num_basis)
        self.barrier = Barrier(explr_space)

        self.lambdak = self.basis.lambdak
        # self.lambdak = np.exp(-0.8*np.linalg.norm(self.basis.k, axis=1))

        self.phik = self.basis.convert_phi2phik(target_dist.grid_vals, target_dist.grid)

        self.u_seq = np.zeros((model.action_space_dim, self.N))
        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), 0.0)))
        # self.u_seq = np.random.uniform(0.0, 0.01, size=(model.action_space_dim,self.N))
        # self.u_def = self.u_seq

        self.R = np.array([[0.01, 0.0],
                           [0.0, 0.001]])

        self.Rinv = np.linalg.inv(self.R)
        self.q = 1.0

        self.pT = np.zeros(model.state_space_dim)

        self.past_states = None


    def controls(self, x):
        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(self.model.action_space_dim)

        # forward simulate trajectory
        xt = np.zeros((self.model.state_space_dim, self.N))
        dfk = []
        fdx = []
        fdu = []
        dbar= []
        for i in range(self.N):
            xt[:,i] = x
            dfk.append(self.basis.dfk(x))
            fdx.append(self.model.fdx(x, self.u_seq[:,i]))
            fdu.append(self.model.fdu(x))
            dbar.append(self.barrier.dx(x))
            x = self.model.step(x, self.u_seq[:,i], self.dt)

        ck = np.zeros(self.basis.tot_num_basis)

        # Average ck based on T not N ???
        # add past states to traj for cks
        if self.past_states is not None:
            past_states = self.past_states.copy()
            total_traj = np.append(past_states, xt, axis=1)
            ck = self.basis.convert_traj2ck(total_traj)
        else:
            ck = self.basis.convert_traj2ck(xt)

        fourier_diff = self.lambdak * (ck - self.phik)


        # backwards pass
        rho = self.pT
        for i in reversed(range(self.N)):
            # Did not add -2*q/N
            edx = np.zeros(self.model.state_space_dim)
            edx[self.model.explr_dim] = self.q * np.dot(fourier_diff, dfk[i])

            bdx = np.zeros(self.model.state_space_dim)
            bdx[self.model.explr_dim] = dbar[i]

            # ergodic metric is not used to update heading
            # rho = rho - self.dt * (-edx - np.dot(fdx[i].T, rho))
            rho = rho - self.dt * (-edx - bdx - np.dot(fdx[i].T, rho))

            # TRY += here ??? (add Udeff in Eqn 8)
            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])


        if self.past_states is None:
            self.past_states = x.reshape(self.model.state_space_dim,1)
        else:
            self.past_states = np.append(self.past_states, x.reshape(self.model.state_space_dim,1), axis=1)

        return self.u_seq[:,0].copy()
