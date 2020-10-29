import numpy as np
import matplotlib.pyplot as plt
import scipy

from fourier import Basis
from target_distribution import TargetDist
from single_integrator import SingleIntegrator
from control import ErgodicControl
from controlSE2 import ErgodicControlSE2
from controlKL import ErgodicControlKL
from barrier import Barrier
from cart import Cart
from omni import Omni

# [xmin = ymin], [xmax = ymax]
explr_space = np.array([0.0, 1.0])
t_dist = TargetDist(num_pts=50)

model = Cart()
# model = Omni()

# erg_ctrl = ErgodicControlSE2(explr_space, model, t_dist, horizon=0.5, num_basis=5)
erg_ctrl = ErgodicControlKL(explr_space, model, t_dist, horizon=0.5, num_samples=10**2)

x_curr = np.array([0.0, 0.0, 0.0])
t_curr = 0.0
tf = 10
dt = 0.1
N = int(tf/dt)
trajectory = np.zeros((3,N))
i = 0


plt.figure(dpi=110,facecolor='w')
xy, vals = t_dist.get_grid_spec()
plt.contourf(*xy, vals, levels=10)

while i < N:
    u = erg_ctrl.controls(x_curr)
    x_curr = model.step(x_curr, u, dt)
    trajectory[:,i] = x_curr
    t_curr  = t_curr + dt
    i = i + 1
    plt.scatter(x_curr[0], x_curr[1])
    plt.pause(0.01)

# plt.figure(dpi=110,facecolor='w')
# xy, vals = t_dist.get_grid_spec()
# plt.contourf(*xy, vals, levels=10)
# plt.scatter(trajectory[0], trajectory[1])
plt.show()
