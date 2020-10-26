import numpy as np
import matplotlib.pyplot as plt
import scipy

from fourier import Basis
from target_distribution import TargetDist
from single_integrator import SingleIntegrator
from control import ErgodicControl
from controlSE2 import ErgodicControlSE2
from controlKL import ErgodicControlKL
from controlMPPIKL import ErgodicControlMPPIKL
from barrier import Barrier
from cart import Cart

# [xmin = ymin], [xmax = ymax]
explr_space = np.array([0.0, 1.0])

#######################################################################
# basis = Basis(explr_space, num_basis=4)

# # normalization factors
# print("hk")
# print(basis.hk.reshape(num_basis,num_basis))
#
# # random trajectory
# xt = np.random.uniform(0.0, 1.0, size=(10,2))
# # print(xt[0])
#
# # basis functions for a point
# print("fk")
# print(basis.fk(xt[0]))
#
# # derivative of basis functions for a point
# print("dfk")
# print(basis.dfk(xt[0]))
#
# print("cks")
# print(basis.convert_traj2ck(xt))


#######################################################################
# basis2 = Basis(explr_space, num_basis=5)
# print(basis2.lambdak)

t_dist = TargetDist(num_pts=50)
# xy, vals = t_dist.get_grid_spec()
# plt.contourf(*xy, vals, levels=10)
# plt.show()
# print(t_dist.grid_vals.shape)
# print(t_dist.grid.shape)

# # phiks for distibution
# phik = basis2.convert_phi2phik(t_dist.grid_vals, t_dist.grid)
# # print(phiks.shape)
# print("phik")
# print(phik)
#
# # Convert phiks back to phi
# phi = basis2.convert_phik2phi(phik, t_dist.grid)
# print("phi")
# print(phi)
#
# print("min phi: ", min(phi))
#
# plt.contourf(phi.reshape(50,50))
# plt.show()


#######################################################################
# model = SingleIntegrator()
# erg_ctrl = ErgodicControl(explr_space, model, t_dist, horizon=0.5, num_basis=5)
#
# x_curr = np.array([0.5, 1.2])
# t_curr = 0.0
# tf = 10
# dt = 0.1
# N = int(tf/dt)
# trajectory = np.zeros((2,N))
# i = 0
#
# # erg_ctrl.controls(t_curr, x_curr)
#
# while i < N:
#     u = erg_ctrl.controls(t_curr, x_curr)
#     x_curr = model.step(x_curr, u, dt)
#     trajectory[:,i] = x_curr
#     t_curr  = t_curr + dt
#     i = i + 1
#
#
# plt.figure(dpi=110,facecolor='w')
# xy, vals = t_dist.get_grid_spec()
# plt.contourf(*xy, vals, levels=10)
# plt.scatter(trajectory[0], trajectory[1])
# plt.show()

# convert trajectory to distibution
# ck = erg_ctrl.basis.convert_traj2ck(trajectory)
# phik_xt = erg_ctrl.basis.convert_ck2dist(ck)
# plt.figure(dpi=110,facecolor='w')
# plt.contourf(*xy, phik_xt.reshape((50,50)), levels=10)
# plt.show()


#######################################################################
# model = Cart()
# # model = SingleIntegrator()
# erg_ctrl = ErgodicControlSE2(explr_space, model, t_dist, horizon=0.5, num_basis=5)
#
# # x_curr = np.array([0.5, 1.2, np.pi/6])
# x_curr = np.array([0.0, 0.0, 0.0])
# t_curr = 0.0
# tf = 20
# dt = 0.1
# N = int(tf/dt)
# trajectory = np.zeros((3,N))
# i = 0
#
# # erg_ctrl.controls(t_curr, x_curr)
#
# plt.figure(dpi=110,facecolor='w')
# xy, vals = t_dist.get_grid_spec()
# plt.contourf(*xy, vals, levels=10)
#
# while i < N:
#     u = erg_ctrl.controls(x_curr)
#     x_curr = model.step(x_curr, u, dt)
#     trajectory[:,i] = x_curr
#     t_curr  = t_curr + dt
#     i = i + 1
#     plt.scatter(x_curr[0], x_curr[1])
#     plt.pause(0.01)
#
# # plt.figure(dpi=110,facecolor='w')
# # xy, vals = t_dist.get_grid_spec()
# # plt.contourf(*xy, vals, levels=10)
# # plt.scatter(trajectory[0], trajectory[1])
# plt.show()

# plt.figure(dpi=110,facecolor='w')
# plt.plot(trajectory[2])
# plt.show()

#######################################################################
model = Cart()
erg_ctrl = ErgodicControlKL(explr_space, model, t_dist, horizon=0.5, num_samples=10**2)

# x_curr = np.array([0.1, 0.9, 0.0])
x_curr = np.array([0.1, 0.0, 0.0])
t_curr = 0.0
tf = 10
dt = 0.1
N = int(tf/dt)
trajectory = np.zeros((model.state_space_dim,N))
i = 0

# erg_ctrl.controls(x_curr)

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

#######################################################################
# model = Cart()
# erg_ctrl = ErgodicControlMPPIKL(explr_space, model, t_dist, horizon=0.5, num_samples=10**2)
#
# # x_curr = np.array([0.5, 0.5, np.pi/6])
# x_curr = np.array([0.1, 0.2, 0.0])
# t_curr = 0.0
# tf = 30
# dt = 0.01
# N = int(tf/dt)
# trajectory = np.zeros((model.state_space_dim,N))
# i = 0
#
# # erg_ctrl.controls(x_curr)
#
# plt.figure(dpi=110,facecolor='w')
# xy, vals = t_dist.get_grid_spec()
# plt.contourf(*xy, vals, levels=10)
#
# while i < N:
#     u = erg_ctrl.controls(x_curr)
#     x_curr = model.step(x_curr, u, dt)
#     trajectory[:,i] = x_curr
#     t_curr  = t_curr + dt
#     i = i + 1
#     plt.scatter(x_curr[0], x_curr[1])
#     plt.pause(0.01)
#
#
# # plt.figure(dpi=110,facecolor='w')
# # xy, vals = t_dist.get_grid_spec()
# # plt.contourf(*xy, vals, levels=10)
# # plt.scatter(trajectory[0], trajectory[1])
# plt.show()

#######################################################################
# grid = np.load('/home/boston/ergodic_exploration_ws/src/ergodic_exploration/data/map.npy')
# plt.figure(dpi=110,facecolor='w')
# plt.plot(grid)
# plt.show()





















#
