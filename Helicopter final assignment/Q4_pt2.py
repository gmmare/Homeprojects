import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import pandas as pd

# ref data
g = 9.81
CL_alpha = 5.7
solidity = 0.075
gamma = 6
C_d_fus = 1.5
m = 2200
rho = 1.225
vtip = 200
r = 7.32
rpm_rad = vtip / r
I_yy = 10615
h = 1
S = np.pi * r ** 2
tau = 0.1

t0 = 0
tmax = 40
dt = 0.1

# initial settings state variables
u_init = 0.0  # velocity X relative to body
w_init = 0.0  # velocity Y relative to body
q_init = 0.0  # pitch rate positive nose up
theta_f_init = 0.00  # fuselage pitch angle positive pitch up
lambda_i_init = np.sqrt(m * 9.81 / (S * 2 * rho)) / vtip  # nondimensionalised instanteneous induced velocity

# initial settings input variables
theta_cyclic_init = 0
theta_collective_init = 6 * np.pi / 180

# setting up state & control vectors
state_init = np.array([u_init,
                       w_init,
                       q_init,
                       theta_f_init,
                       lambda_i_init])

control_init = np.array([theta_cyclic_init,
                         theta_collective_init])

# Creating time range for maneuver
t_range = np.arange(t0, tmax + dt, dt)

# initializing arrays for storing state data
state_range = np.array([state_init])
control_range = np.array([control_init])
lam_i = []
a1_check = []
phi_check = []

for i, t in enumerate(t_range[:-1]):

    state_v = state_range[-1]
    control_v = control_range[-1]

    # getting state & control variables
    u = state_v[0]
    w = state_v[1]
    q = state_v[2]
    theta_f = state_v[3]
    lambda_i = state_v[4]

    #attitude control:
    theta_f_wish = 0

    #cyclic laws
    if t >= 0.5 and t < 1:
        theta_cyclic = 1 * np.pi / 180

    else:
        theta_cyclic = 0

    if t >= 15.0:
        theta_cyclic = 0.2 * (theta_f - theta_f_wish) + 0.2 * q

    #Collective laws
    theta_collective = 6 * np.pi / 180

    # checking quadrant angle
    if u == 0:
        if w > 0:
            phi = np.pi / 2
        else:
            phi = -np.pi / 2
    else:
        phi = np.arctan(w/u)

        if u < 0:
            phi = phi + np.pi

    V = np.sqrt((u ** 2) + (w ** 2))
    alpha_c = theta_cyclic - phi

    mu = (V / vtip) * np.cos(alpha_c)
    lambda_c = (V / vtip) * np.sin(alpha_c)
    qdiml = q/rpm_rad

    teller = -16/gamma * qdiml + 8/3*mu*theta_collective - 2*mu*(lambda_c+lambda_i)
    a1 = teller/(1-.5*mu**2)

    # CT calculation #correct working
    CT_bem = 0.25 * CL_alpha * solidity * ((2 / 3) * theta_collective * (1 + (3 / 2) * mu ** 2) -
                                           (lambda_c + lambda_i))

    CT_glau = 2 * lambda_i * np.sqrt(((V / vtip) * np.cos(alpha_c - a1)) ** 2 +
                                 ((V / vtip) * np.sin(alpha_c - a1) + lambda_i) ** 2)

    T = CT_bem * rho * vtip ** 2 * S

    #eq of motion
    u_dot = -g * np.sin(theta_f) - 0.5 * (C_d_fus / m) * rho * u * V + \
            (T / m) * np.sin(theta_cyclic - a1) - q * w

    w_dot = g * np.cos(theta_f) - 0.5 * (C_d_fus / m) * rho * w * V - \
            (T / m) * np.cos(theta_cyclic - a1) + q * u

    # g*cos(pitch(i))-cds/mass*.5*rho*w(i)*vv(i)-thrust(i)/mass*cos(helling(i))+q(i)*u(i)

    q_dot = - (T / I_yy) * h * np.sin(theta_cyclic - a1)

    theta_f_dot = q

    # updating state vector
    state_v_new = np.array([u + u_dot * dt,
                            w + w_dot * dt,
                            q + q_dot * dt,
                            theta_f + theta_f_dot * dt,
                            lambda_i + ((CT_bem - CT_glau) / tau) * dt])



    # updating state vector
    control_v_new = np.array([theta_cyclic,
                              theta_collective])

    state_range = np.append(state_range, [state_v_new], axis=0)
    control_range = np.append(control_range, [control_v_new], axis=0)

    lam_i.append(lam_i)
    a1_check.append(a1)
    phi_check.append(phi)


# plt.plot(t_range[:-1], a1_check)
plt.plot(t_range, state_range[:, 0], label='u')
plt.plot(t_range, state_range[:, 3] * 180 / np.pi, label='theta_f')
# plt.plot(t_range, control_range[:,1], label='longitudinal')
# plt.plot(t_range, control_range[:,0], label='cycli')
plt.legend()
plt.show()
