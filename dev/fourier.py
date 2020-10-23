import numpy as np

class Basis(object):
    def __init__(self, explr_space, num_basis=5):
        self.L = explr_space[1] - explr_space[0]
        self.n = 2
        self.tot_num_basis = num_basis**self.n

        k = np.meshgrid(*[[i for i in range(num_basis)] for _ in range(self.n)])
        self.k = np.c_[k[0].ravel(), k[1].ravel()]

        # normalization factors
        self.hk = np.zeros(self.tot_num_basis)
        for i, k in enumerate(self.k):
            if np.prod(k) < 1e-5:
                self.hk[i] = 1.
            else:
                top = np.prod(self.L * (2.0 * k * np.pi + np.sin(2.0 * k *np.pi)))
                bot = 16.0 * np.prod(k) * np.pi**2
                self.hk[i] = top/bot

        # lambdaks
        self.lambdak = np.zeros(self.tot_num_basis)
        for i, k in enumerate(self.k):
            self.lambdak[i] = (1.0 / (1.0 + self.k[i,0]**2 + self.k[i,1]**2))**((self.n+1.0)/2.0)


    def fk(self, x):
        """ x is a point (x, y)
        Output: is length K^n

        Note: normalization by 1/hk is removed
        """
        fkvec = np.zeros(self.tot_num_basis)
        for i, k in enumerate(self.k):
            fkvec[i] = np.cos(k[0] * np.pi/self.L * x[0]) * np.cos(k[1] * np.pi/self.L * x[1])
        return fkvec
        # return np.prod(np.cos(np.pi *x / self.L * self.k), 1)


    def dfk(self, x):
        """ x is a point (x, y)
        Output: is length 2 [DF/dx DF/dy]

        Note: normalization by 1/hk is removed
        """
        dfkvec = np.zeros((self.tot_num_basis, self.n))
        for i, k in enumerate(self.k):
            k1term = k[0] * np.pi/self.L
            k2term = k[1] * np.pi/self.L

            DFkx1 = -k1term * np.sin(k1term * x[0]) * np.cos(k2term * x[1])
            DFkx2 = -k2term * np.cos(k1term * x[0]) * np.sin(k2term * x[1])

            dfkvec[i,:] = np.array([DFkx1, DFkx2])

        return dfkvec


    def convert_traj2ck(self, xt):
        """
        Converts a trajectory into its time-averaged
        statistics in the Fourier domain

        Note: N is the number of the steps in the trajectory and is used
              instead of T the length in time as noted in the papers

        Note: assumes xt is a row vector wrt to forward simulation
        """
        N = xt.shape[1]
        # print("convert_traj2ck (lenght of xt): ", N)
        mat = np.zeros((self.tot_num_basis, N))
        for i in range(N):
            mat[:,i] = self.fk(xt[:,i])

        # sum across columns
        return (1.0/N) * np.sum(mat, axis=1)


    def convert_phi2phik(self, phi_val, phi_grid):
        """
        Converts the distribution to the fourier decompositions
        """
        num_pts = len(phi_val)
        mat = np.zeros((self.tot_num_basis, num_pts))
        for i in range(num_pts):
            mat[:,i] = self.fk(phi_grid[i,:]) * phi_val[i]

        # sum across columns
        return np.sum(mat, axis=1)
        # return np.sum([self.fk(x) * v for v, x in zip(phi_val, phi_grid)], axis=0)


    def convert_phik2phi(self, phik, phi_grid):
        """
        Reconstructs phi from the Fourier terms
        """
        num_pts = phi_grid.shape[0]
        phi = np.zeros(num_pts)
        for i in range(num_pts):
            phi[i] = np.dot(self.fk(phi_grid[i,:]), phik)

        # phi_val = np.stack([np.dot(basis.fk(x), phik) for x in phi_grid])
        return phi


    def convert_ck2dist(self, ck, grid=None, size=1.0):
        """
        Converts a ck into its time-averaged
        statistics
        """
        grid = np.meshgrid(*[np.linspace(0, size)
                                for _ in range(self.n)])
        grid = np.c_[grid[0].ravel(), grid[1].ravel()]

        num_pts = grid.shape[0]
        phi = np.zeros(num_pts)
        for i in range(num_pts):
            phi[i] = np.dot(self.fk(grid[i,:]), ck)
        return phi









#
