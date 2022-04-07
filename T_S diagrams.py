import matplotlib.pyplot as plt
import math

cpa = 1000
cpg = 1150
R = 287

#temps
T2 = 288.15
T21 = 333.66
T25 = 379.2
T3 = 834.6
T4 = 1627
T45 = 1235.3
T5 = 937.8
T8 = 804.99

#pressures
p2 = 99298.5
p21 = 158877.6
p3 = 3098113.2
p4 = 2974188.672
p45 = 875674.1455
p5 = 257719.3442
p8 = 137377.0492
s_tab = [0]
t_tab = [T2, T21, T3, T4, T45, T5 ,T8]
p_tab = [p2, p21, p3, p4, p45, p5, p8]
def entropy_a(Tt1, Tt2, P1, P2):
    ds = cpa * math.log(Tt2/Tt1, math.e) - R * math.log(P2/P1, math.e)

    return ds

def entropy_g(Tt1, Tt2, P1, P2):
    ds = cpg * math.log(Tt2/Tt1, math.e) - R * math.log(P2/P1, math.e)

    return ds

for i in range(2):
    Tt1 = t_tab[i]
    Tt2 = t_tab[i + 1]
    P1 = p_tab[i]
    P2 = p_tab[i + 1]
    print(entropy_a(Tt1, Tt2, P1, P2) + s_tab[-1])
    s_tab.append(entropy_a(Tt1, Tt2, P1, P2) + s_tab[-1])

for i in range(2, 6, 1):
    Tt1 = t_tab[i]
    Tt2 = t_tab[i + 1]
    P1 = p_tab[i]
    P2 = p_tab[i + 1]
    print(entropy_g(Tt1, Tt2, P1, P2) + s_tab[-1])
    s_tab.append(entropy_a(Tt1, Tt2, P1, P2) + s_tab[-1])

plt.plot(s_tab, t_tab, marker=".")

plt.ylabel('T [K]')
plt.xlabel('S [J/K]')
plt.title('T-S diagram of CFM56-5B at sea level static conditions')
plt.show()
