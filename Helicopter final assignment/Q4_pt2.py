import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import pandas as pd

#parameters
m = 3600    #kg
n = 4       #blades
r = 5.5     #m
rpm = 475   #rpm
max_speed = 278 #kmh
C_d_fus = 0.15 # [-]
rho = 1.225
g = 9.81
S = 2 * 1.7
c = 0.35
h = 1.5
rpm = 375
rpm_rad = rpm * 2 * np.pi / 60
I_yy = 25000
solidity = n * c / (np.pi * r) #sigma
CL_alpha = 2*np.pi
gamma = (rho * CL_alpha * c * r**2)/I_yy

# Built space state like model where the control input depends on distance to target attitude velocity
# state vector = [u, w, q, theta_f, lambda_i], input vector = [theta_cyclic, theta_collective]
def GetThrustCoeff(state_v, control_v):

    u = state_v[0]
    w = state_v[1]
    q = state_v[2]
    lambda_i = state_v[4]
    theta_cyclic = control_v[0]
    theta_collective = control_v[1]

    #calculating constants
    V = np.sqrt(u**2 + w**2)
    alpha_c = theta_cyclic - np.arctan2(w, u)
    mu = V/(rpm_rad * r) * np.cos(alpha_c)
    lambda_c = V/(rpm_rad * r) * np.sin(alpha_c)

    a1 = ((8/3) * mu * theta_collective - 2 * mu * (lambda_c + lambda_i) -(16/gamma) * (q/rpm_rad))

    #calculating CT
    CT_bem = (1/4) * CL_alpha * solidity * ((2/3) * theta_collective * (1 + (3/2) * mu**2)
                                            - (lambda_c + lambda_i))

    CT_glau = lambda_i * np.sqrt((V/(rpm_rad * r) * np.cos(alpha_c - a1))**2 +
                                 (V/(rpm_rad * r) * np.sin(alpha_c - a1) + lambda_i)**2)

    return CT_bem, CT_glau, V, a1

#simulation settings
t0 = 0
tmax = 50
dt = 0.1

#initial settings state variables
u_init = 0                  #velocity X relative to body
w_init = 0                  #velocity Y relative to body
q_init = 0                  #pitch rate positive nose up
theta_f_init = 0            #fuselage pitch angle positive pitch up
lambda_i_init = np.sqrt(m*9.81/(S*2*rho))/(rpm_rad*r)  #nondimensionalised instanteneous induced velocity

#initial settings input variables
theta_cyclic_init = 0
theta_collective_init = 0

#setting up vectors
state_init = np.array([u_init,
                       w_init,
                       q_init,
                       theta_f_init,
                       lambda_i_init])
control_init = np.array([theta_cyclic_init,
                         theta_collective_init])

# Creating time range for maneuver
t_range = np.arange(t0, tmax + dt, dt)

#initializing arrays for storing state data
state_range = np.array([])
control_range = np.array([])
 
for i, t in enumerate(t_range[:-1]):
    if t == 0:
        state_v = state_init
        control_v = control_init
    else:
        state_v = state_range[-1]
        control_v = control_range[-1]

    # change in inst. induced velocity
    tau = 0.1
    CT_bem, CT_glau, V, a1 = GetThrustCoeff(state_v, control_v)
    dvi_dt = (CT_bem - CT_glau) / tau

    #getting state & control variables
    u = state_v[0]
    w = state_v[1]
    q = state_v[2]
    theta_f = state_v[3]
    lambda_i = state_v[4]
    theta_cyclic = control_v[0]
    theta_collective = control_v[1]

    #calculating thrust & drag
    T = dvi_dt* rho * (rpm_rad * r)**2 * S
    D = C_d_fus * (1/2) * rho * V**2 * S

    #eq of motion
    u_dot = -g * np.sin(theta_f) - (D/m) * (u/V) + (T/m) * np.sin(theta_cyclic - a1) - q * w
    w_dot = g * np.cos(theta_f) - (D/m) * (w/V) + (T/m) * np.cos(theta_cyclic - a1) - q * u
    q_dot = -(T/I_yy) * h * np.sin(theta_cyclic - a1)
    theta_f_dot = q

    #updating state vector
    state_v_new = np.array([u + u_dot * dt,
                            w + w_dot * dt,
                            q + q_dot * dt,
                            theta_f + theta_f_dot * dt,
                            lambda_i + dvi_dt * dt])

    #Control laws PID part



