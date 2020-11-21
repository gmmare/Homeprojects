import matplotlib.pyplot as plt
import math

import numpy as np





# ref data
f=open('F1_10.dat', 'r')
data =[]
for line in f.readlines():
    data.append(line)
f.close()

# sim data cf
sim_cf=open('cpdat.dat', 'r')
data_cf =[]
for line in sim_cf.readlines():
    data_cf.append(line)
sim_cf.close()

# sim data cf
sim_cp=open('cpdat.dat', 'r')
data_cp =[]
for line in sim_cp.readlines():
    data_cp.append(line)
sim_cp.close()

#======================= reference data reading =======================
x_data_cp = []
cp_data = []

# succ side cp
for i in range(35, 72):
    # x position of the data
    num_x = float(data[i][14:20])
    exp_x = float(data[i][21:24])
    x_data_cp.append(round((num_x * (10 ** exp_x)),4))

    # cp data
    num_cp= float(data[i][27:33])
    exp_cp = float(data[i][35:38])
    cp_data.append(round((num_cp * (10 ** exp_cp)),4))


# cf side
x_data_cf = []
cf_data = []

for i in range(74, 87):
    # x position of the data
    num_x = float(data[i][14:20])
    exp_x = float(data[i][21:24])
    x_data_cf.append(round((num_x * (10 ** exp_x)),4))

    # cp data
    num_cf= float(data[i][28:33])
    exp_cf = float(data[i][35:38])
    cf_data.append(round((num_cf * (10 ** exp_cf)),4))

# ======================= sim data reading =======================
x_sim_cp = []
cp_sim = []


# succ side cp sim
for i in range(62, 105):

    data_line = data_cp[i]
    data_line.strip(' ')
    data_split = data_line.split(',')

    # x position of the data
    num_x = float(data_split[0][0:10])
    exp_x = float(data_split[0][11:14])
    x_sim_cp.append(round((num_x * (10 ** exp_x)),4))

    # cp data
    start_split = data_split[1].strip('\n')
    # print(data_split[1], 'after strip')
    cp = start_split.split('e')

    num_cp = float(cp[0])
    exp_cp = float(cp[1])
    cp_sim.append(round((num_cp * (10 ** exp_cp)),4))


# # cf side
# x_sim_cf = []
# cf_sim = []
#
# for i in range(74, 87):
#     # x position of the data
#     num_x = float(data_cf[i][14:20])
#     exp_x = float(data_cf[i][21:24])
#     x_sim_cf.append(round((num_x * (10 ** exp_x)),4))
#
#     # cp data
#     num_cf= float(data_cf[i][28:33])
#     exp_cf = float(data_cf[i][35:38])
#     cf_sim.append(round((num_cf * (10 ** exp_cf)),4))

#plotting Cp
plt.plot(x_data_cp, cp_data, label='ref')
plt.plot(x_sim_cp, cp_sim, label='sim')
plt.xlabel('Chordwise position (x/c) [-]')  # Add an x-label to the axes.
plt.ylabel('Cp [-]')  # Add a y-label to the axes.
plt.legend()  # Add a legend.
plt.show()

# # #plotting Cf
# plt.plot(x_data_cf, cf_data, label='test1')
# plt.plot(x_data_cf, cf_data, label='test2')
# plt.xlabel('Chordwise position (x/c) [-]')  # Add an x-label to the axes.
# plt.ylabel('Cf [-]')  # Add a y-label to the axes.
# plt.legend()  # Add a legend.
# plt.show()

