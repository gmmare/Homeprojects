import matplotlib.pyplot as plt
import math

def tw(P_ratio):
    mflow = 196.3714878
    Pamb = 101325
    Fan_isen = 0.85
    nozzle_eff = 0.98
    k_a = 1.4
    cp_a = 1000
    v0 = 102.0522219
    T_in = 293.184
    Pin = 107719.9653
    m_pl = 2225
    m_train = 70
    motor_eff = 0.95
    pmu_eff = 0.9
    bat_eff = 0.85
    inverter_eff = 0.95

    T_out = T_in*(1 + (1/Fan_isen) * ((P_ratio**((k_a -1)/k_a)) -1))

    t7_t8 = T_out*nozzle_eff*(1 - (Pamb/(P_ratio*Pin))**((k_a-1)/k_a))

    vjet = math.sqrt(2 * cp_a*t7_t8)

    T_gen = mflow*(vjet - v0)

    v_ind = T_gen / (2 * mflow)

    P = T_gen * v_ind                     #actuator disk method

    E_fan = 600 * P /(motor_eff * bat_eff * pmu_eff * inverter_eff)

    m_bat = E_fan / (3600 * 200)

    t_w = T_gen/((m_train + m_bat + (m_pl/4))*9.81)
    return t_w, m_bat

xtab = []
ytab = []
mbat = []
dpr = 0.01
for i in range(15):
    xtab.append(1 + i * dpr)
    ytab.append(tw(1 + i * dpr)[0])
    mbat.append(tw(1 + i * dpr)[1])


plt.plot(xtab, mbat)
plt.vlines(1.078,0,mbat[-1])
plt.ylabel('battery mass [kg]')
plt.xlabel('Fan pressure ratio [-]')
plt.show()

plt.hlines(1.1,xtab[0],xtab[-1], label='T/W = 1.1')
plt.xlabel('Fan pressure ratio [-]')
plt.ylabel('T/W ratio [-]')
plt.plot(xtab, ytab)
plt.legend()
plt.show()