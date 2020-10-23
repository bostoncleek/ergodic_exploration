import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, trapz

# Domain of world
X1min = 0.0           # lower bound on x1 (x-axis)
X1max = 1.0            # upper bund on x1 (x-axis)
X2min = 0.0           # lower bound on x2 (y-axis)
X2max = 1.0            # upper bund on x2 (y-axis)

L1 = X1max - X1min      # length of axis x1
L2 = X2max - X2min      # length of axis x2
L = L1 = L2

n = 2                   # dims in ergodic exploration
K = 5                   # basis functions along each axis
# Nk = K**n               # total number of basis functions
res = 0.1               # discretization along each axis

# Length of trajectory
T = 5.0
t0 = 0.0
N = 100
tspan = [t0, T]
dt = float(T/(N-1))
tvec = np.linspace(min(tspan), max(tspan), N)

# Ergodic metric cost
q = 1.0


Q = np.array([[1.0, 0.0],
               [0.0,  1.0]])

R = np.array([[1.0, 0.0],
               [0.0,  1.0]])


u0 = np.array([-0.5, -0.5])
x0 = np.array([1.0, 0])
z0 = np.array([0.0, 0.0])


def spatial_distribution(x):
    """ Multivariate normal distribution
    Args:
    x: (x1, x2)
    """
    mu = np.array([0.5, 0.5])
    sigma = np.eye(2)
    temp = np.dot(np.dot(x - mu, np.linalg.inv(sigma)), x - mu)
    return (np.linalg.det(2.0 * np.pi * sigma)**(-0.5)) * np.exp(-0.5 * temp)


def lambdak(k1, k2):
    """ Compose constant """
    s = float((n+1)/2.0)
    return (1.0 + (k1**2 + k2**2))**(s)


def hk(k1, k2):
    """  Fourier basis normalization factor
    Args:
    k1: basis number
    k2: basis number
    """

    if (k1 * k2) < 1e-5:
        return 1.0

    else:
        x1term = L1 * (2.0 * k1 * np.pi + np.sin(2.0 * k1 * np.pi))
        x2term = L2 * (2.0 * k2 * np.pi + np.sin(2.0 * k2 * np.pi))

        return (x1term * x2term) / (16 * k1 * k2 * np.pi**2)

    # return np.sqrt(L1*L2)


def Fk(x, h, k1, k2):
    """ Fourier basis function of agents's state
    Args:
    x: state (x, y) correspond to (x1, x2)
    h: fourier basis normalization factor
    k1: basis number
    k2: basis number
    """
    # return (1.0/h) * np.cos(k1 * np.pi/L1 * x[0]) * np.cos(k2 * np.pi/L2 * x[1])
    return np.cos(k1 * np.pi/L1 * x[0]) * np.cos(k2 * np.pi/L2 * x[1])


def DFk(x, h, k1, k2):
    """ Derivative of fourier basis
    Args:
    x: state(x, y) correspond to (x1, x2)
    h: fourier basis normalization factor
    k1: basis number
    k2: basis number
    """
    k1term = k1 * np.pi/L1
    k2term = k2 * np.pi/L2

    DFkx1 = -k1term * np.sin(k1term * x[0]) * np.cos(k2term * x[1])
    DFkx2 = -k2term * np.cos(k1term * x[0]) * np.sin(k2term * x[1])

    # return (1.0/h) * np.array([DFkx1, DFkx2])
    return np.array([DFkx1, DFkx2])


def Ck(x, h, k1, k2):
    """ Fourier coeffecient of spatial statics of trajectory
    Args:
    x: trajectory (x, y)
    h: normalization
    k1: basis number
    k2: basis number
    """
    coeff = np.zeros(N)
    for i in range(N):
        coeff[i] = Fk(x[:,i], h, k1, k2)
    return (1.0/N) * trapz(coeff, tvec)


def phik(h, k1, k2):
    """ Fourier coefficient of the spatial distribution
    Args:
    h_k: normalization factor
    k1: basis number
    k2: basis number
    """
    # Discretize domain
    l1_domain = np.arange(X1min, X1max, res)
    l2_domain = np.arange(X2min, X2max, res)

    coeff = 0.0
    for i in l1_domain:
        for j in l2_domain:
            pt = np.array([i, j])
            coeff = coeff + spatial_distribution(pt) * Fk(pt, h, k1, k2)
    return coeff * res * res


def hkMat():
    """ 2D matrix of fourier bases normalization factors """
    hk_mat = np.zeros((K,K))
    for k1 in range(K):
        for k2 in range(K):
            hk_mat[k1, k2] = hk(k1, k2)
    return hk_mat

#     k = np.meshgrid(*[[i for i in range(K)] for _ in range(n)])
#     k = np.c_[k[0].ravel(), k[1].ravel()]

#     hk = np.zeros(k.shape[0])

#     for i, k in enumerate(k):
#         if np.prod(k) < 1e-5:
#             hk[i] = 1.
#         else:
#             top = np.prod(L * (2.0 * k * np.pi + np.sin(2.0 * k *np.pi)) )
#             bot = 16.0 * np.prod(k) * np.pi**2
#             hk[i] = top/bot

#     return hk.reshape(K,K)


def phikMat(hk_mat):
    """ 2D matrix of spatial coefficients
    Args:
    hk_mat: normalization factors
    """
    phik_mat = np.zeros((K,K))
    for k1 in range(K):
        for k2 in range(K):
            h = hk_mat[k1, k2]
            phik_mat[k1, k2] =  phik(h, k1, k2)
    return phik_mat


#################################################
# Pre compute hk and phik
hk_mat = hkMat()
phik_mat = phikMat(hk_mat)
#################################################

print("hks: \n", hk_mat)
print("phiks: \n", phik_mat)


def ergodic_metric(x):
    """ How far a trajectory is from being ergodic with respect to a distribution
    Args:
    x: trajectory (x, y)
    """
    erg = 0.0
    for k1 in range(K):
        for k2 in range(K):
            # normalization
            h = hk_mat[k1, k2]

            # trajectory fourier coefficient
            c_k = Ck(x, h, k1, k2)

            # distribution fourier coefficient
            phi_k = phik_mat[k1, k2]

            erg = erg + lambdak(k1, k2) * (c_k - phi_k)**2
    return erg


def get_a(x, index):
    """ Compose a(t).T
    Args:
    x: trajectory (x, y)
    index: position in trajectory used to compose DFk
    """
    a = np.zeros((1,2))
    for k1 in range(K):
        for k2 in range(K):
            # normalization
            h = hk_mat[k1, k2]

            # trajectory fourier coefficient
            c_k = Ck(x, h, k1, k2)

            # distribution fourier coefficient
            phi_k = phik_mat[k1, k2]

            a = a + (c_k - phi_k) * DFk(x[:,index], h, k1, k2)

    return (((2.0 * q) / N) * a).flatten()


def xdot(t, x, u):
    """ Kinematic model xdot = [u1, u2]
    Args:
    t: current time
    x: (x, y)
    u: (u1, u2)
    """
    i = int(round((N-1)*(t-t0)/(T-t0)))
    if (i == N or i > N):
        print("WARNING i is out of bounds")
        i = N-1

    return np.array([u[0,i], u[1,i]])


def getA():
    """ State linearization A = D1f(x,u) """
    return np.zeros((2,2))


def getB():
    """ Control linearization B = D2f(x,u) """
    return np.eye(2)


def zdot(t, z, v):
    """ zdot = Az + Bv
    Args:
    t: current time
    z: state perturbations
    v: control perturbations
    """
    i = int(round((N-1)*(t-t0)/(T-t0)))
    if (i == N or i > N):
        print("WARNING i is out of bounds")
        i = N-1

    A = getA()
    B = getB()
    return np.dot(A,z) + np.dot(B,v[:,i])


def J(x, u):
    """ Cost function
    Args:
    x: trajectory (x, y, theta)
    u: constrol signal (u1, u2)
    """
    # Ergodic cost
    erg_cost = q * ergodic_metric(x)

    # Control cost
    integrand = np.zeros(N)
    for i in range(N):
        integrand[i] = np.dot(np.dot(u[:,i], R), u[:,i])

    return erg_cost + 0.5 * trapz(integrand, tvec)


def DuJ(x, u, z, v):
    """ Direction derivative cost wrt to controls
    Args:
    x: trajectory (x, y)
    u: constrol signal (u1, u2)
    z: trajectory perturbations
    v: control perturbations
    """
    integrand = np.zeros(N)
    for i in range(N):
        a = get_a(x, i)
        b = np.dot(u[:,i], R)
        integrand[i] = np.dot(a, z[:,i]) +  np.dot(b, v[:,i])
    return trapz(integrand, tvec)


def line_search(x, u, z, v, alpha=0.4, beta=0.1):
    """ Performs Armijo line search
    Args:
    x: trajectory (x, y)
    u: constrol signal (u1, u2)
    z: trajectory perturbations
    v: control perturbations
    """

    i = 0;
    converged = False
    cost = J(x, u)
    DJ_zeta = DuJ(x, u, z, v)

    while not converged:
        gamma = beta**i
        i = i + 1

        # peturb controls
        unew = u + gamma*v

        # perturb trajectory
        solx = solve_ivp(xdot, tspan, x0, t_eval=tvec, args=(unew,))
        xnew = solx.y

        new_cost = J(xnew, unew)
        threshold = cost +  alpha*gamma*DJ_zeta

#         print("Must be less than: ", threshold)
#         print("Cost: ", new_cost)

        if (new_cost < threshold or np.isclose(new_cost, threshold)):
            converged = True

    print("Armjo Converged in", i, "steps")
    return xnew, unew


def zdot2(t, z, x, u, Pvec, rvec):
    """ Compose xdot given solution to riccati equations
    Args:
    t: current time
    z: current z in integration
    x: trajectory (x, y)
    u: constrol signal (u1, u2)
    Pvec: P
    rvec: r
    """

    i = int(round((N-1)*(t-t0)/(T-t0)))
    if (i == N or i > N):
        print("WARNING i is out of bounds")
        i = N-1

    A = getA()
    B = getB()

    Rinv = np.linalg.inv(R)
    b = np.dot(R, u[:,i])

    # index into P and r
    j = int((N-1) - i)
    P = Pvec[:,j].reshape(n,n)
    r = rvec[:,j]

    v = -Rinv.dot(B.T).dot(P).dot(z) - Rinv.dot(B.T).dot(r) - Rinv.dot(b)
    return A.dot(z) + B.dot(v)


def get_v(x, u, zvec, Pvec, rvec):
    """ Compose v given z
    Args:
    x: trajectory (x, y)
    u: constrol signal (u1, u2)
    zvec: trajectory perturbations
    Pvec: P
    rvec: r
    """

    v = np.zeros((2,N))
    for i in range(N):
        A = getA()
        B = getB()

        Rinv = np.linalg.inv(R)
        b = np.dot(R, u[:,i])
        z = zvec[:,i]

        # index into P and r
        j = int((N-1) - i)
        P = Pvec[:,j].reshape(n,n)
        r = rvec[:,j]

        v[:,i] = -Rinv.dot(B.T).dot(P).dot(z) - Rinv.dot(B.T).dot(r) - Rinv.dot(b)
    return v


def riccati(t, inputs, x, u):
    """ Solves Riccati equations for P and r
    Args:

    """

    P = inputs[0:4].reshape(n,n)
    r = inputs[4:].reshape(n)

    i = int(round((N-1)*(t-t0)/(T-t0)))
    if (i == N or i > N):
        print("WARNING i is out of bounds")
        i = N-1

    A = getA()
    B = getB()

    Rinv = np.linalg.inv(R)
    a = get_a(x, i)
    b = np.dot(R, u[:,i])

    Pdot = -(P.dot(A) + A.T.dot(P) - P.dot(B).dot(Rinv).dot(B.T).dot(P) + Q)
    rdot = -((A - B.dot(Rinv).dot(B.T).dot(P)).T.dot(r) + a - P.dot(B).dot(Rinv).dot(b))

    return np.concatenate((Pdot.flatten(), rdot))


def descent_direction(x, u):
    """ Compose z and v
    Args:
    x: trajectory (x, y)
    u: constrol signal (u1, u2)
    """

    # 1) solve Riccati from T -> t0
    # P(T) = 0
    PT = np.zeros((n,n))
    #r(T) = 0
    rT = np.zeros(n)
    # Terminal condition for riccati eqns
    inputT = np.concatenate((PT.flatten(), rT.flatten()))


    solPr = solve_ivp(riccati, [T, t0], inputT, t_eval=tvec[::-1], args=(x, u))
    Pvec = solPr.y[0:4,:]
    rvec = solPr.y[4:,:]

    # 2) solve zdot from t0 -> T
    solz = solve_ivp(zdot2, [t0, T], z0, t_eval=tvec, args=(x, u, Pvec, rvec))
    z = solz.y

    # 3) solve for v using z
    v = get_v(x, u, z, Pvec, rvec)

    return z, v


def iLQR(x, u, eps):
    """ Ergodic exploration using iLQR
    Args:
    x: initial trajectory
    u: initial control signal
    eps: stopping criteria
    """
    erg = []
    cost = []
    erg.append(ergodic_metric(x))
    cost.append(J(x, u))

    i = 0
    tol_reached = False
    while not tol_reached:
        print(i)

        # 1) Descent direction
        z, v = descent_direction(x, u)

        # 2) Armijo line search
        x, u = line_search(x, u, z, v)

        # 3) Perturbation norm
#         norm = np.linalg.norm(np.vstack((z, v)))
        norm = np.linalg.norm(DuJ(x, u, z, v))
        print("norm of zeta: ", np.linalg.norm(np.vstack((z, v))))
        print("norm of DJzeta: ", norm)

        erg.append(ergodic_metric(x))
        cost.append(J(x, u))

        i = i + 1

        if (i == 20):
            break

        if (norm < eps):
            tol_reached = True

    return x, u, cost, erg




# initial controls
init_uvec = np.vstack((np.full((1,N), u0[0]), np.full((1,N), u0[1])))
# init_uvec = np.vstack((np.full((1,N), -0.02), 0.5*np.sin(tvec) ))

# initial trajectory
solx = solve_ivp(xdot, tspan, x0, t_eval=tvec, args=(init_uvec,))
init_xvec = solx.y


xvec, uvec, costvec, ergvec = iLQR(init_xvec, init_uvec, 0.01)


plt.figure(dpi=110,facecolor='w')
plt.plot(init_xvec[0], init_xvec[1])
plt.plot(xvec[0], xvec[1])
plt.legend(("Initial Trajectory", "Optimal Trajectory"))
plt.xlabel(r'x(t)')
plt.ylabel(r'y(t)')
plt.grid(True)
plt.show()

plt.figure(dpi=110,facecolor='w')
plt.plot(costvec)
plt.ylabel("Cost")
plt.xlabel("Iteritaion")
plt.grid(True)
plt.show()

plt.figure(dpi=110,facecolor='w')
plt.plot(ergvec)
plt.ylabel("Ergodic Metric")
plt.xlabel("Iteritaion")
plt.grid(True)
plt.show()


# plt.figure(dpi=110,facecolor='w')
# plt.plot(tvec, uvec[0])
# plt.plot(tvec, uvec[1])
# plt.legend([r'$u_{1}(t)$', r'$u_{2}(t)$'])
# plt.xlabel(r't')
# plt.ylabel(r'u(t)')
# plt.grid(True)
# plt.show()




# p1 = np.array([0.0, 0.0])
# p2 = np.array([1.0, 0.0])
# p3 = np.array([0.0, 1.0])
# p4 = np.array([5.0, 5.0])
#
# print(spatial_distribution(p1))
# print(spatial_distribution(p2))
# print(spatial_distribution(p3))
# print(spatial_distribution(p4))

# p1 = np.array([0.0, 0.0])
# phi_val = spatial_distribution(p1)
#
# phi_grid = np.meshgrid(*[np.linspace(0, 1.0, int(np.sqrt(len(phi_val))))
#                         for _ in range(2)])
#
# phi_grid = np.c_[phi_grid[0].ravel(), phi_grid[1].ravel()]
#
# # assert phi_grid.shape[0] == phi_val.shape[0], 'samples are not the same'
# #
# # return np.sum([basis.fk(x) * v for v, x in zip(phi_val, phi_grid)], axis=0)
#
# print(phi_grid)
