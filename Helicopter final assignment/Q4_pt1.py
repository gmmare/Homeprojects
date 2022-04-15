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
C_dS = 1.2 # [-]
rho = 1.225
S = 2 * 1.7
c = 0.35
rpm = 375
rpm_rad = rpm * 2 * np.pi / 60

#trim calculations
def Get_Trim(V):
    W = m * 9.81
    D = C_dS * 1/2 * rho * V**2

    #calculating thrust coeff
    T = np.sqrt(W**2 + D**2)
    C_t = T/(rho * (rpm_rad * r) ** 2 * np.pi * r ** 2)
    mu = V/(rpm_rad * r)

    #solving for lambda_i (non dimensionalised induced velocity)
    x = sy.symbols('lambda_i')
    ans = sy.solve(C_t - 2 * x * sy.sqrt((mu * np.cos(D / W)) ** 2
                                        + (mu * np.sin(D / W)
                                           + x) ** 2))
    lambda_i = ans[0]


    #calculating other variables

    solidity = n * c / (np.pi * r)
    CL_alpha = 2*np.pi

    #setting up matrix for solving
    A11 = 1 + (3/2) * mu**2
    A12 = -(8/3) * mu
    A21 = -mu
    A22 = (2/3) + mu**2

    B1 = -2 * mu**2 * D/W - 2 * mu * lambda_i
    B2 = (4/solidity) * (C_t/CL_alpha) + mu * (D/W) + lambda_i

    A = [[A11, A12],
         [A21, A22]]
    B = [[B1], [B2]]

    x_trim = np.linalg.inv(A).dot(B)
    x_trim = x_trim * 180/np.pi

    return x_trim[0][0], x_trim[1][0]

trim_state = pd.DataFrame(columns=['V','theta_c','theta_o'])
for V in range(100):
    t_c, t_o = Get_Trim(V)
    trim_state.loc[V] = [V, t_c, t_o]


plt.plot(trim_state['V'], trim_state['theta_c'], label='Cyclic pitch')
plt.plot(trim_state['V'], trim_state['theta_o'], label='collective pitch')
plt.xlabel("Forward velocity [m/s]")
plt.ylabel("Control actuation [degr]")
plt.legend()
plt.show()
print("done")
