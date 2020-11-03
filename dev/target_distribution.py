import numpy as np

class TargetDist(object):
    def __init__(self, num_pts=50):
        self.num_pts = num_pts

        # create 2D grid on domain [0 1] x [0 1]
        grid = np.meshgrid(*[np.linspace(0, 1, num_pts) for _ in range(2)])
        self.grid = np.c_[grid[0].ravel(), grid[1].ravel()]

        # means and variance of each distirbution

        # 3 distributions
        # self.means = [np.array([0.7, 0.7]), np.array([0.3,0.3]), np.array([0.2,0.8])]
        # self.vars  = [np.array([0.1,0.1])**2, np.array([0.2,0.2])**2, np.array([0.4,0.4])**2]

        # 2 distributions
        self.means = [np.array([0.7, 0.7]), np.array([0.3,0.3])]
        self.vars  = [np.array([0.1,0.1])**2, np.array([0.1,0.1])**2]

        # 1 distribution
        # self.means = [np.array([0.7, 0.7])]
        # self.vars = [np.array([0.1, 0.1])**2]

        self.grid_vals = self.__call__(self.grid)


    def sample_points(self, x):
        return self.__call__(x)

    def sample_grid_spec(self, x, sample_vals):
        xy = []
        num_pts = int(np.sqrt(sample_vals.shape))
        for g in x.T:
            xy.append(
                np.reshape(g, newshape=(num_pts, num_pts))
            )
        return xy, sample_vals.reshape(num_pts, num_pts)


    def get_grid_spec(self):
        xy = []
        for g in self.grid.T:
            xy.append(
                np.reshape(g, newshape=(self.num_pts, self.num_pts))
            )
        return xy, self.grid_vals.reshape(self.num_pts, self.num_pts)


    def __call__(self, x):
        val = np.zeros(x.shape[0])
        for m, v in zip(self.means, self.vars):
            innerds = np.sum((x-m)**2 / v, 1)
            val += np.exp(-innerds/2.0)

        # normalizes the distribution
        val /= np.sum(val)
        return val













#
