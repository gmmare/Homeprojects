import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import pandas as pd

#parameters
M = 3600    #kg
n = 4       #blades
r = 5.5     #m
rpm = 475   #rpm
max_speed = 278 #kmh
C_d_fus = 0.15 # [-]
rho = 1.225
S = 2 * 1.7
c = 0.35
rpm = 375

###============ Q2 induced velocity plot ============###
def get_vi(V_bar): #approximation


    # lowsspeed flight:
    coeffs = [1, V_bar**2, -1]
    x = np.roots(coeffs)

    v_i = np.sqrt(x[1])

    return v_i

def get_vi2(V_forward): #approximation
    W = M * 9.81
    D = C_d_fus * 1/2 * rho * S * V_forward**2

    #calculating alpha
    alpha_d = D/W
    T = W/np.cos(alpha_d)

    #non dimensionalising
    V_ih = np.sqrt(T/(2 * rho * np.pi * r**2))
    V_bar = V_forward/V_ih

    #calculating velocity components
    vcos = V_bar * np.cos(alpha_d)
    vsin = V_bar * np.sin(alpha_d)

    #solving num, first defining the variable that will be solved for
    v_i_bar = sy.symbols('v_i_bar')
    ans = sy.solve(v_i_bar**2 * (vcos**2 + (vsin - v_i_bar)**2) - 1)

    return V_bar, -ans[0]

#getting numerical points
x_range = []
y_range = []
for i in range(0, 60):
    x, y = get_vi2((i))
    x_range.append(x)
    y_range.append(y)

#getting points using approximations
v_bar_range = np.arange(0,5,0.1)
v_i_range = []

for i, V in enumerate(v_bar_range):
    v_i_range.append(get_vi(V))

#plotting
plt.scatter(v_bar_range, v_i_range, color='r', linewidth=0, label='Approximation')
plt.plot(x_range, y_range, label='Solved numerically')
ax = plt.gca()
plt.xlabel(r'$\bar{V}$ [-] Forward flight')
plt.ylabel(r'$v_i$ [-] Induced velocity')
plt.legend()
plt.show()

###============ Q3 Performance calculation ============###
W = M * 9.81
P_i = W * np.sqrt(W/(2 * rho * np.pi * r**2))
print("Ideal power P_i=", P_i * 1e-3, "[kW]")

#BEM power in hover
r_effective = 0.97 * r
solidity = n * c / (np.pi * r)
tip_speed = rpm * (2 * np.pi/60) * r_effective
C_d_profile = 0.012
k_factor = 1.1
extra_loss = 1.1


P_hov_bem = (P_i * k_factor + (C_d_profile/8) * rho * solidity * (tip_speed ** 3) * np.pi * (r_effective**2)) * extra_loss

#ACT power in hover
FM = 0.725
P_hov_act = P_i/FM * extra_loss

print("P_hov with ACT = ", P_hov_act * 1e-3, "[kW]")
print("P_hov with BEM = ", P_hov_bem * 1e-3, "[kW]")

#BEM power forward flight

def get_power(V_forward):
    to_kw = 1e-3
    tail_rotor = 1.1
    tip_speed_ratio = V_forward/tip_speed
    P_par = S * C_d_fus * 0.5 * rho * V_forward ** 3 * to_kw
    P_profdrag = (solidity * C_d_profile / 8) * rho * (tip_speed**3) * np.pi * (r_effective ** 2) \
    * (1 + 4.65 * tip_speed_ratio**2) * to_kw
    V_bar, v_i = get_vi2(V_forward)
    P_ind = k_factor * W * v_i * np.sqrt(W/(2 * np.pi * rho * r**2)) * to_kw

    return P_par, P_profdrag, P_ind, (P_par + (P_profdrag + P_ind)*tail_rotor)

P_res = pd.DataFrame(columns=['V','Ppar','Ppd','Pind','Ptot'])

for V in range(0, 100):
    P_par, P_pd, P_ind, P_tot = get_power(V)
    P_res.loc[V] = [V, P_par, P_pd, P_ind, P_tot]

#converting data frame to all floats
P_res = P_res.astype(float)

#plotting
plt.plot(P_res['V'], P_res['Ppar'], label=r'$P_{parasitic}$')
plt.plot(P_res['V'], P_res['Ppd'], label=r'$P_{profile - drag}$')
plt.plot(P_res['V'], P_res['Pind'], label=r'$P_{induced}$')
plt.plot(P_res['V'], P_res['Ptot'], label=r'$P_{total}$')
plt.xlabel('Forward velocity [m/s]')
plt.ylabel('Power [kW]')
plt.legend()
plt.show()

#max endurance
V_endurance = P_res['Ptot'].idxmin()
print('Endurance speed:', V_endurance, '[m/s], Ptot =', P_res['Ptot'][V_endurance], '[kW]')

#max range
def Get_Performance(V_range, P_range, V_endurance, P_endurance):
    V_range = V_range.to_numpy()
    P_range = P_range.to_numpy()
    tangent_diff = []

    for i, V in enumerate(V_range[1:-2]):
        dx_local = V_range[i+1] - V_range[i]
        dy_local = P_range[i+1] - P_range[i]

        tangent_local = dy_local/dx_local

        dx_origin = V
        dy_origin = P_range[i]

        tangent_origin = dy_origin/dx_origin
        tangent_diff.append(abs((tangent_local-tangent_origin)/tangent_origin))

    #getting the tangent point
    min_value = min(tangent_diff)
    min_index = tangent_diff.index(min_value)
    print("Speed for maximum range:", V_range[min_index], "[m/s], Ptot =", P_range[min_index], '[kW]')

    #plotting
    plt.plot(V_range, P_range)
    plt.plot([0, V_range[min_index], V_range[-1]], [0, P_range[min_index], V_range[-1] * (P_range[min_index]/V_range[min_index])])
    plt.plot([V_endurance, V_endurance], [0, P_endurance], '--')
    plt.xlabel('Forward velocity [m/s]')
    plt.ylabel('Power [kW]')
    plt.show()

    return

Get_Performance(P_res['V'], P_res['Ptot'], V_endurance, P_res['Ptot'][V_endurance])








