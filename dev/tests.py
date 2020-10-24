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

# plt.figure(dpi=110,facecolor='w')
# plt.plot(trajectory[2])
# plt.show()

#######################################################################
model = Cart()
erg_ctrl = ErgodicControlKL(explr_space, model, t_dist, horizon=0.5, num_samples=100**2)

# x_curr = np.array([0.5, 1.2, np.pi/6])
x_curr = np.array([0.5, 0.5, 0.0])
t_curr = 0.0
tf = 20
dt = 0.1
N = int(tf/dt)
trajectory = np.zeros((3,N))
i = 0

# erg_ctrl.controls(x_curr)

while i < N:
    u = erg_ctrl.controls(x_curr)
    x_curr = model.step(x_curr, u, dt)
    trajectory[:,i] = x_curr
    t_curr  = t_curr + dt
    i = i + 1


plt.figure(dpi=110,facecolor='w')
xy, vals = t_dist.get_grid_spec()
plt.contourf(*xy, vals, levels=10)
plt.scatter(trajectory[0], trajectory[1])
plt.show()

# plt.figure(dpi=110,facecolor='w')
# plt.plot(trajectory[2])
# plt.show()


# u = np.zeros((10, 2))
# u[0,0] = 0.5
# u[0,1] = 0.7
# u[9,0] = 1
# u[9,1] = 2
#
# print(u)
# print("---------")
# print(u[:-1,:])
# print("---------")
# print(u[1:,:])
# print("---------")
# u[:-1,:] = u[1:,:]
# print(u)
#
# print("---------")
# print(u[-1,:])
# print("---------")
# print(u[-1])
# u[-1,:] = u[-1]
# print("---------")
# print(u)




#######################################################################
# model = Cart()
# x_eq = np.array([0.0, 0.0, 0.707])
# # x_eq = np.array([0.0, 0.0, np.pi/4.0])
# u_eq = np.array([0.0, 0.0])
# #
# A = model.fdx(x_eq, u_eq)
# B = model.fdu(x_eq)
#
# A = np.zeros((3,3))
# B = np.eye(3)
#
# print("Dynamics better be 0: ", model.f(x_eq, u_eq))
# print(A)
# print(B)
#
# Q = np.eye(3)
# R = np.eye(3)
#
# P = scipy.linalg.solve_discrete_are(A, B, Q, R)
# K = np.linalg.inv(B.T.dot(P).dot(B) + R).dot(B.T.dot(P).dot(A))
#
# print("P ", np.round(P,3))
# print("Optimal gain: ", np.round(K,3))





#
