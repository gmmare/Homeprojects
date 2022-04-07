import matplotlib.pyplot as plt
import numpy as np
import sympy as sy

# plt.rcParams['text.usetex'] = True

#parameters
M = 3600    #kg
n = 4       #blades
r = 5.5     #m
rpm = 475   #rpm
max_speed = 278 #kmh
C_d_fus = 0.12 # [-]
rho = 1.225
S = 2 * 1.7



# ###============ Q2 induced velocity plot ============###
# def get_vi(V_bar): #approximation
#
#
#     # lowsspeed flight:
#     coeffs = [1, V_bar**2, -1]
#     x = np.roots(coeffs)
#
#     v_i = np.sqrt(x[1])
#
#     return v_i
#
# def get_vi2(V_forward): #approximation
#     W = M * 9.81
#     D = C_d_fus * 1/2 * rho * S * V_forward**2
#
#     #calculating alpha
#     alpha_d = D/W
#     T = W/np.cos(alpha_d)
#
#     #non dimensionalising
#     V_ih = np.sqrt(T/(2 * rho * np.pi * r**2))
#     V_bar = V_forward/V_ih
#
#     #calculating velocity components
#     vcos = V_bar * np.cos(alpha_d)
#     vsin = V_bar * np.sin(alpha_d)
#
#     #solving num, first defining the variable that will be solved for
#     v_i_bar = sy.symbols('v_i_bar')
#     ans = sy.solve(v_i_bar**2 * (vcos**2 + (vsin - v_i_bar)**2) - 1)
#
#     return V_bar, -ans[0]
#
# #getting numerical points
# x_range = []
# y_range = []
# for i in range(0, 60):
#     x, y = get_vi2((i))
#     x_range.append(x)
#     y_range.append(y)
#
# #getting points using approximations
# v_bar_range = np.arange(0,5,0.1)
# v_i_range = []
#
# for i, V in enumerate(v_bar_range):
#     v_i_range.append(get_vi(V))
#
# #plotting
# plt.scatter(v_bar_range, v_i_range, color='r', linewidth=0, label='Approximation')
# plt.plot(x_range, y_range, label='Solved numerically')
# ax = plt.gca()
# plt.xlabel(r'$\bar{V}$ [-] Forward flight')
# plt.ylabel(r'$v_i$ [-] Induced velocity')
# plt.legend()
# plt.show()

###============ Q3 Performance calculation ============###
W = M * 9.81
P_i = W * np.sqrt(W/(2 * rho * np.pi * r**2))
print("ideal power P_i:", P_i * 1e-3, "kW")


