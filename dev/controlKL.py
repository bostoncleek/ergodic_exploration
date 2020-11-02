import numpy as np
import matplotlib.pyplot as plt

from barrier import Barrier

class ErgodicControlKL(object):
    def __init__(self, explr_space, model, target_dist, horizon, num_samples):
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
        self.R = np.array([[0.5, 0.0],
                           [0.0, 0.01]])

        self.Rinv = np.linalg.inv(self.R)

        # Penalty on ergodic measure
        self.q = 1.0

        # Uncertainty in (x,y) position
        self.sigma = np.eye(2) * 0.1
        self.sigmaInv = np.linalg.inv(self.sigma)

        # self.eta = np.linalg.det(2.0*np.pi*self.sigma)**0.5

        # Init condition for co-state integration p(T) = 0
        self.pT = np.zeros(model.state_space_dim)


    # def ergodic_measure_terms(self, s, ps, x):
    #     # pdot based on eqn 20
    #     sum = np.zeros(self.model.state_space_dim)
    #     for i in range(self.num_samples):
    #         sum[self.model.explr_dim]  +=  2.0 * ps[i] * np.dot(s[:,i] - x[self.model.explr_dim], self.sigmaInv)
    #     # Not sure about the negaive sign but seems to drive it in the right direction ????
    #     return -1.*sum


    def ergodic_measure_terms2(self, si, xt):
        """
        Compose q and dq using a reimann sum
        Args:
            si: point (x,y) uniformly sampled from exploration domain
            xt: trajectory
        Returns:
            q: time average trajectory statistics (Eqn 3)
            dq: derivative of q w.r.t. the state
        """

        q = 0.0
        dq = np.zeros(self.model.state_space_dim)
        for i in range(self.N):

            # NOTE: self.model.explr_dim = [0, 1], the exploration domain
            # is less than the dimension of the robot's state space.
            # The robot states space inludes heading but the ergodic measure does not
            # consider the heading.
            q_integrand = np.exp(-0.5 * np.dot(si - xt[self.model.explr_dim, i], self.sigmaInv).dot(si - xt[self.model.explr_dim, i]))
            dq_integrand = q_integrand * np.dot(si - xt[self.model.explr_dim, i], self.sigmaInv)

            q += q_integrand
            dq[self.model.explr_dim] += dq_integrand

        # q = (1./(self.N * self.eta)) * q
        # dq = (1./(self.N * self.eta)) * dq

        return q, -1.*dq


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

        # KL-E^3 Base Algo. Line 11
        # sample points in (x,y) space
        # these are used to approximate the KL divergence
        # KL divergence approximates the ergodic measure
        s = np.random.uniform(self.explr_space[0], self.explr_space[1], (2, self.num_samples))
        ps = self.target_dist.sample_points(s.T)

        # xy, vals = self.target_dist.sample_grid_spec(s.T, ps)
        # plt.contourf(*xy, vals, levels=10)
        # plt.show()

        # edx is the derivative of the ergodic measure approximated by KL-div.
        # w.r.t. to the state
        # This is the first terms in Eqn 10.
        edx = np.zeros(self.model.state_space_dim)
        for i in range(self.num_samples):
            q, dq = self.ergodic_measure_terms2(s[:,i], xt)
            edx += (ps[i]/q) * dq

        # Apply penaliation to ergodic measure
        edx = self.q * edx
        # print(edx)

        # KL-E^3 Base Algo. Lines 12-18
        # Solve for the co-state variable by integrating backwards in time
        # backwards pass
        rho = self.pT
        for i in reversed(range(self.N)):
            # Collect the boundary derivative
            bdx = np.zeros(self.model.state_space_dim)
            bdx[self.model.explr_dim] = dbar[i]

            ## SIGN ERROR  on edx ?????
            # Update the co-state variable
            rho = rho - self.dt * (-edx -bdx -np.dot(fdx[i].T, rho))

            # Update the control sequence
            self.u_seq[:,i] = -np.dot(np.dot(self.Rinv, fdu[i].T), rho)

            # Normalize controls if there are too large
            if (np.abs(self.u_seq[:,i]) > 1.0).any():
                self.u_seq[:,i] /= np.linalg.norm(self.u_seq[:,i])

        # KL-E^3 Base Algo. Line 19
        # Send the first controls to actuators
        # print(self.u_seq)
        return self.u_seq[:,0].copy()
