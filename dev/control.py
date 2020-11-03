import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from fourier import Basis
from barrier import Barrier

class ErgodicControl(object):
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


        # self.u_seq = np.zeros((2,self.N))
        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), 0.0)))
        self.u_seq = np.random.uniform(0.0, 0.01, size=(2,self.N))
        # self.u_def = self.u_seq

        self.R = np.array([[1.0, 0.0],
                           [0.0, 1.0]])

        self.Rinv = np.linalg.inv(self.R)
        self.q = 2.0

        self.pT = np.zeros(2)

        self.past_states = None


    # def pdot(self, t, p, fourier_diff, dfk, fdx, t0, tf):
    #     i = int(round((self.N-1)*(t-t0)/(tf-t0)))
    #     # print("time: ", round(t,3), " index: ", i)
    #     p = -np.dot(fourier_diff, dfk[i]) - np.dot(fdx[i].T, p)
    #     return p


    def controls(self, t_curr, x_curr):
        x = x_curr

        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(2)

        # forward simulate trajectory
        xt = np.zeros((2, self.N))
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
            edx = self.q * np.dot(fourier_diff, dfk[i])
            # edx = (2.0 * self.q / self.N) * np.dot(fourier_diff, dfk[i])

            # rho = rho - self.dt * (-edx - np.dot(fdx[i].T, rho))
            rho = rho - self.dt * (-edx - dbar[i] - np.dot(fdx[i].T, rho))

            # TRY += here ??? (add Udeff in Eqn 8)
            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])


        if self.past_states is None:
            self.past_states = x.reshape(2,1)
        else:
            self.past_states = np.append(self.past_states, x.reshape(2,1), axis=1)

        return self.u_seq[:,0].copy()



    # def controls(self, t_curr, x_curr):
    #     # 1) forward simulate trajectory
    #     t0 = t_curr
    #     tf = t_curr + self.horizon
    #
    #     tvec = np.linspace(t0, tf, self.N)
    #     xt = solve_ivp(self.model.f, [t0, tf], x_curr, t_eval=tvec, args=(self.u_seq, t0, tf, self.N)).y
    #
    #     # some required terms
    #     dfk = []
    #     fdx = []
    #     fdu = []
    #     for i in range(self.N):
    #         dfk.append(self.basis.dfk(xt[:,i]))
    #         fdx.append(self.model.fdx())
    #         fdu.append(self.model.fdu())
    #
    #
    #     # 2) fourier coeffecients of trajectory
    #     ck = self.basis.convert_traj2ck(xt)
    #     # print("ck shape: ",ck.shape)
    #
    #     # difference in fourier coeffecients
    #     fourier_diff = self.lambdak * (ck - self.phik)
    #     # print("fourier_diff")
    #     # print(fourier_diff)
    #
    #     # for i in range(self.N):
    #     #     print(np.dot(fourier_diff, dfk[i]))
    #
    #
    #     # 3) p(t) from tf -> t0
    #     revtec = tvec[::-1]
    #     p = solve_ivp(self.pdot, [tf, t0], self.pT, t_eval=revtec, args=(fourier_diff, dfk, fdx, t0, tf)).y
    #
    #     # print("p")
    #     # print(p)
    #     # plt.figure(dpi=110,facecolor='w')
    #     # plt.plot(p[0])
    #     # plt.plot(p[1])
    #     # # plt.plot(self.u_seq[0])
    #     # # plt.plot(self.u_seq[1])
    #     # plt.grid(True)
    #     # plt.show()
    #
    #
    #     # 4) optimal u
    #     for i in range(self.N):
    #         j = int((self.N-1) - i)
    #         self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), p[:,j])
    #         # print(np.dot(np.dot(self.Rinv, fdu[i].T), p[:,j]))
    #     # print(self.u_seq)
    #
    #     # send first controls to robot
    #     u_apply = self.u_seq[:,0]
    #
    #     # shift controls over
    #     self.u_seq[:,:-1] = self.u_seq[:,1:]
    #
    #     # set last controls to zero
    #     self.u_seq[:,-1]  = np.zeros(2)
    #
    #     return u_apply






























#
