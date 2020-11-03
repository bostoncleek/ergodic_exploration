import numpy as np
import matplotlib.pyplot as plt
import scipy

from target_distribution import TargetDist
from controlSE2 import ErgodicControlSE2
from controlKL import ErgodicControlKL
from controlKLIVP import ErgodicControlKLIVP
from barrier import Barrier
from cart import Cart
from omni import Omni

# [xmin = ymin], [xmax = ymax]
explr_space = np.array([0.0, 1.0])

# PDF representing a 2D distribution
# where the robot will find the most usesful information
t_dist = TargetDist(num_pts=50)

# 2D Kinematic cart
# model = Cart()
model = Omni()

# This produces the desired output
# erg_ctrl = ErgodicControlSE2(explr_space, model, t_dist, horizon=0.5, num_basis=5)

# This does not work as expected
erg_ctrl = ErgodicControlKL(explr_space, model, t_dist, horizon=0.5, num_samples=10**2, buffer_size=50)

# erg_ctrl = ErgodicControlKLIVP(explr_space, model, t_dist, horizon=0.5, num_samples=10**2, buffer_size=50)


x_curr = np.array([0.5, 0.5, 0.0])
t_curr = 0.0
tf = 50
dt = 0.1
N = int(tf/dt)
trajectory = np.zeros((3,N))
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


#########################################
# x_eq = np.zeros(3)
# u_eq = np.zeros(4)

# x_eq = np.array([0.0, 1.0, np.pi/2.0])
# u_eq = np.array([1.0, -1.0, 1.0, -1.0])
#
# omni = Omni()
# A = omni.fdx(x_eq, u_eq)
# B = omni.fdu(x_eq)
# print("f: ", omni.f(x_eq, u_eq))
# print(A)
# print(B)
#
# Q = np.eye(3)
# R = np.eye(4)*0.1
#
# P = scipy.linalg.solve_discrete_are(A, B, Q, R)
# print(P)
#
# K = np.linalg.inv(B.T.dot(P).dot(B) + R).dot(B.T.dot(P).dot(A))
# print(K)



#
