import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

Q = np.array([[1e3, 0.0, 0.0],
              [0.0, 1e3, 0.0],
              [0.0, 0.0, 1.0]])

R = np.array([[0.1, 0.0],
               [0.0,  0.1]])

P1 = np.array([[1e3, 0.0, 0.0],
               [0.0, 1e3, 0.0],
               [0.0, 0.0, 1e3]])


# # wheel radius
# radius = 0.033
# # wheel base
# wheel_base = 0.16
# # max motor rotational speed
# u_max = 6.0



def euler(f, x, u, dt):
    """ Euler integration for one time step"""
    x = np.copy(x)
    u = np.copy(u)
    new = x + f(x, u)*dt
    return new


def cart(x, u):
    """ Diff Drive Robot
    Args:
        x: (x,y,theta)
        u: u[0] = ul (left wheel velocity), u[1] = ur (right wheel velocity)
    """
    # res = np.array([(radius/2.0)*np.cos(x[2])*(u[0] + u[1]),
    #                 (radius/2.0)*np.sin(x[2])*(u[0] + u[1]),
    #                 (radius/wheel_base)*(u[1] - u[0])])
    # return res
    return np.array([np.cos(x[2])*u[0], np.sin(x[2])*u[0], u[1]])


# def ref_traj(t):
# #     return np.array([radius*(u_max/2.0)*t, 0.0, np.pi/2.0])
#     return np.array([(2.0/np.pi)*t, 0.0, np.pi/2.0])



def loss(x, xd, u):
    """ Compose loss
    Args:
        x: (x,y,theta)
        xd: desired pose (xd,yd,thetad)
        u: u[0] = ul (left wheel velocity), u[1] = ur (right wheel velocity)
    """
    return 0.5*((x-xd).T.dot(Q).dot(x-xd) + u.T.dot(R).dot(u))


def d_goal(x):
    """Compose distance to waypoint
    Args:
        x: (x,y,theta)
    """
    return np.linalg.norm(xT[0:2]-x[0:2])


def MPPI(x0, u0, dt, cov, N, K, lam, dthresh, t_curr):
    """
    Model predictive path integral control for kinematic cart waypoint following
    Args:
        x0: initial pose (x0,y0,theta0)
        u0: initial control signal for left and righ wheel velocities (shape = (2 x N))
        dt: time difference
        cov: covariance matrix used for sampling perturbations (shape = (2x2))
        N: number of time steps used for forward propagation of model in time horizon
        K: number of rollouts
        lam: temperature paramter scales cost
        dthresh: distance threshold to waypoint
        t_curr: current time
    """
    xinit = np.copy(x0)
    uvec = np.copy(u0)

    # control signal and trajectory
    u_opt = []
    xvec = []

    # store perturbuations
    du1 = np.zeros((N,K))
    du2 = np.zeros((N,K))
    # record cost at each discretization for each rollout
    J = np.zeros((N,K))

    converged = False
    cnt = 0
    while not converged:
        for k in range(K):
            # IC for each rollout
            x = xinit
            # perturbuations
#             du = np.random.multivariate_normal([0.0, 0.0], cov, N).T
            duL = np.random.normal(0.0, cov[0,0], N)
            duR = np.random.normal(0.0, cov[1,1], N)
            du = np.vstack((duL, duR))

            du1[:,k] = du[0]
            du2[:,k] = du[1]

            for i in range(N):
                # perturb controls
                u_pert = uvec[:,i] + du[:,i]

                # propagate kinematics
                x = euler(cart, x, u_pert, dt)

                # compose loss
#                 xd = ref_traj(t_curr)
                xd = xT
                J[i,k] = loss(x, xd, u_pert)

            # terminal cost
            J[N-1,k] += (x-xT).T.dot(P1).dot(x-xT)


        # sum across each rollout backwards in time
        # consider how current actions effect future states
        J = np.cumsum(J[::-1], 0)[::-1,:]
#         print(J)

        for i in range(N):
            # min across all rollouts for this point in time
            beta = np.amin(J[i])
            J[i] -= beta

            # compose weights
            w = np.exp(-J[i]/lam) + 1e-8
            w /= np.sum(w)

            # update controls
            uvec[0,i] += np.sum(np.dot(w, du1[i,:]))
            uvec[1,i] += np.sum(np.dot(w, du2[i,:]))


#         print(J)
        # print(uvec)

        # saturate controls
#         uvec = np.clip(uvec, -u_max, u_max)

        # filter controls
        uvec = savgol_filter(uvec, N-1, 3, axis=1)

        # send first control to actuators
        xinit = euler(cart, xinit, uvec[:,0], dt)
        xvec.append(xinit.tolist())
        u_opt.append(uvec[:,0].tolist())

        # update controls for next iteration
        uvec[:,0:N-1] = uvec[:,1:N]
        # update last control
        uvec[:,N-1] = u0[:,N-1]


        d = d_goal(xinit)
        if (cnt % 100 == 0):
            print("iteration:", cnt, " t: ", t_curr, " d: ", d)

        if (d < dthresh or np.isclose(d,dthresh)):
            converged = True

        t_curr += dt
        cnt += 1


    xvec = np.array(xvec).T
    u_opt = np.array(u_opt).T
    return xvec, u_opt



T = 1.0 # horizon
N = 10 # number of time steps
dt = float(T/float(N))

cov = np.array([[0.1, 0.0],
                 [0.0, 0.1]])


# u0 = np.vstack((np.full((1,N), 3.0), np.full((1,N), 0.8)))
u0 = np.zeros((2,N))

x0 = np.array([0.0, 0.0, 0.0])

xT = np.array([1.0, 0.0, np.pi/2.0])



xvec, u_opt = MPPI(x0, u0, dt, cov, N, K=5, lam=0.01, dthresh=0.1, t_curr=0.0)
print("Final pose:", xvec[:,-1])


plt.figure(dpi=110,facecolor='w')
plt.plot(xvec[0], xvec[1])
plt.xlabel(r'x(t)')
# plt.xlim([-0.05, 0.05])
plt.ylabel(r'y(t)')
plt.grid(True)
plt.show()


plt.figure(dpi=110,facecolor='w')
plt.plot(u_opt[0])
plt.plot(u_opt[1])
plt.xlabel(r'iteration')
plt.ylabel(r'u(t)')
plt.legend([r'$u_{l}(t)$', r'$u_{r}(t)$'])
plt.grid(True)
plt.show()
