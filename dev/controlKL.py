import numpy as np
import matplotlib.pyplot as plt

from barrier import Barrier

class ErgodicControlKL(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples, buffer_size):
        # Exporation space domain [0, 1]
        self.explr_space = explr_space
        # Distribution used to compose KL divergence
        self.target_dist = target_dist
        # Control time horizon
        self.horizon = horizon
        # Time step in integration
        self.dt = 0.1
        # Number of pointes in integration
        self.N = int(self.horizon / self.dt)
        # Number of points uniformly sampled to approx. KL divergence
        self.num_samples = num_samples
        # number of past states to store
        self.buffer_size = buffer_size


        # Kinematic model of robot
        self.model = model
        # Barrier Keeps robot from leaving search domain
        self.barrier = Barrier(explr_space)

        # Conrol signal
        self.u_seq = np.zeros((model.action_space_dim, self.N))
        # self.u_seq = np.vstack((np.full((1,self.N), 0.1), np.full((1,self.N), 0.0)))

        # Penalty on controls
        # R[0,0] penalizes forward velocity
        # R[1,1] penalizes rotational velocity
        # self.R = np.array([[0.1, 0.0],
        #                    [0.0, 0.01]])
        # self.R = np.eye(2) * 0.001
        self.R = np.eye(4) * 0.001


        self.Rinv = np.linalg.inv(self.R)

        # Penalty on ergodic measure
        # self.q = 1.0

        # Uncertainty in (x,y) position
        self.sigma = np.eye(2) * 0.0025
        self.sigmaInv = np.linalg.inv(self.sigma)

        # Init condition for co-state integration p(T) = 0
        self.pT = np.zeros(model.state_space_dim)

        self.past_states = None
        self.stored = 0


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
        """
        Compose new control
        Args:
            x: current robot state

        Return:
            u: new control (case: cart this is the linear and angular velocities)
        """

        # shift controls over
        self.u_seq[:,:-1] = self.u_seq[:,1:]

        # set last controls to zero
        self.u_seq[:,-1]  = np.zeros(self.model.action_space_dim)

        # KL-E^3 Base Algo. Lines 3-10
        # forward simulate trajectory
        xt = np.zeros((self.model.state_space_dim, self.N))
        fdx = []
        fdu = []
        dbar= []
        for i in range(self.N):
            xt[:,i] = x
            # jacobain of the Kinematics w.r.t to the state
            fdx.append(self.model.fdx(x, self.u_seq[:,i]))
            # jacobian of the Kinematics w.r.t to the controls
            fdu.append(self.model.fdu(x))
            # barrier derivative
            dbar.append(self.barrier.dx(x))
            x = self.model.step(x, self.u_seq[:,i], self.dt)

        total_traj = None
        if self.past_states is not None:
            total_traj = np.append(self.past_states.copy(), xt, axis=1)
        else:
            total_traj = xt

        # KL-E^3 Base Algo. Line 11
        # sample points in (x,y) space
        # these are used to approximate the KL divergence
        # KL divergence approximates the ergodic measure
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        # xy, vals = self.target_dist.sample_grid_spec(s.T, ps)
        # plt.contourf(*xy, vals, levels=10)
        # plt.show()

        q = np.zeros(self.num_samples)
        for j in range(self.num_samples):
            q[j] = self.traj_stat(s[:,j], total_traj)
        q /= np.sum(q)

        # KL-E^3 Base Algo. Lines 12-18
        # Solve for the co-state variable by integrating backwards in time
        # backwards pass
        rho = self.pT
        for i in reversed(range(self.N)):
            # edx is the derivative of the ergodic measure approximated by KL-div.
            # w.r.t. to the state
            # This is the first terms in Eqn 10.
            edx = np.zeros(self.model.state_space_dim)
            for j in range(self.num_samples):
                # q = self.traj_stat(s[:,j], total_traj)
                dq = self.traj_stat_deriv(s[:,j], xt[:,i])
                edx += (ps[j]/q[j]) * dq
                # edx += (ps[j]/q) * dq

            # edx *= self.q
            # edx[self.model.explr_dim] += 1e-5
            # print(edx)


            # Collect the boundary derivative
            bdx = np.zeros(self.model.state_space_dim)
            bdx[self.model.explr_dim] = dbar[i]

            # Update the co-state variable
            rho = rho - self.dt * (edx -bdx -np.dot(fdx[i].T, rho))

            # Update the control sequence
            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            # Normalize controls if there are too large
            if (np.abs(self.u_seq[:,i]) > 1.0).any() or (np.abs(self.u_seq[:,i]) < -1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

            # print(self.u_seq[:,i])

        if self.past_states is None:
            self.past_states = x.reshape(self.model.state_space_dim,1)
        elif self.stored < self.buffer_size :
            self.past_states = np.append(self.past_states, x.reshape(self.model.state_space_dim,1), axis=1)
            self.stored += 1
        else:
            self.past_states[:,:-1] = self.past_states[:,1:]
            self.past_states[:,-1] = x


        # KL-E^3 Base Algo. Line 19
        # Send the first controls to actuators
        # print(self.u_seq)
        return self.u_seq[:,0].copy()
