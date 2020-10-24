import numpy as np

class TargetDist(object):
    def __init__(self, num_pts=50):
        self.num_pts = num_pts

        # create 2D grid on domain [0 1] x [0 1]
        grid = np.meshgrid(*[np.linspace(0, 1, num_pts) for _ in range(2)])
        self.grid = np.c_[grid[0].ravel(), grid[1].ravel()]
        # two distributions
        # means and variance of each
        self.means = [np.array([0.7, 0.7]), np.array([0.3,0.3])]
        self.vars  = [np.array([0.1,0.1])**2, np.array([0.1,0.1])**2]

        # self.means = [np.array([0.5, 0.5])]
        # self.vars = [np.array([0.1, 0.1])]

        self.grid_vals = self.__call__(self.grid)


    # def spatial_distribution(self, pt):
    #     mu = np.array([0.5, 0.5])
    #     sigma = np.eye(2)
    #     temp = np.dot(np.dot(pt - mu, np.linalg.inv(sigma)), pt - mu)
    #     return (np.linalg.det(2.0 * np.pi * sigma)**(-0.5)) * np.exp(-0.5 * temp)
    #
    # def __call__(self, x):
    #     val = np.zeros(x.shape[0])
    #     for i, p in enumerate(x):
    #         val[i] = self.spatial_distribution(p)
    #     return val

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
