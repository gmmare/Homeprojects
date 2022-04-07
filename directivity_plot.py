import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from scipy import signal

#defining a function to get the average of the last 5% of data points

#converting power spectral density to spl
def pa_to_db1(pa):
    spl = 10 * np.log10(pa / ((2 * 10 ** -5) ** 2)) + 10 * np.log10(800)
    return spl

def pa_to_db2(pa):
    spl = 10 * np.log10(pa / ((2 * 10 ** -5) ** 2)) + 10 * np.log10(1600)
    return spl

def pa_to_db3(pa):
    spl = 10 * np.log10(pa / ((2 * 10 ** -5) ** 2)) + 10 * np.log10(3200)
    return spl

#converting spl to ospl for a pressure band
def compute_ospl(db_list):
    integral = []
    for i in range(len(db_list)):
        integral.append(10 ** (db_list[i]/10))

    total_pa = 10 * np.log10(sum(integral))

    return total_pa

# function for plotting pressure time series
def plot_series(filename):
    with open(filename) as f:
        lines = f.readlines()
    n_last = len(lines)
    y_pa = []
    t = []

    #extracting time series data
    for i in range(5, n_last-3):
        data = lines[i].split(' ')
        if len(data)<2:
            data = [0, y_pa[-1]]
        y_pa.append(float(data[1]))
        t.append(float(data[0]))

    plt.plot(t, y_pa)
    plt.title("theta 2 pi")
    plt.show()

    return

# plot_series("mic_data2/FWH_test.fwh_sim0.txt")

#function to perform fourier transform, split into frequency bands and convert to db
def fourier_filter(filename):
    with open(filename) as f:
        lines = f.readlines()

    #preparing data for fourier transofrm
    n_last = len(lines)
    y_pa = []

    #extracting time series data
    for i in range(5, n_last):
        data = lines[i].split(' ')
        if len(data)<2:
            data = [0, y_pa[-1]]
        y_pa.append(float(data[1]))

    mean = sum(y_pa)/len(y_pa)
    y_pa = [x - mean for x in y_pa]
    dt = 0.2 / len(y_pa)

    # use welch approach for fourier transform
    freq_range, spectrum = signal.welch(y_pa, fs=1/dt)
    spectrum = spectrum * 2

    # prepping data for filtering using pandas dataframe and converting Pa to db
    data_set = pd.DataFrame(spectrum, index=freq_range, columns=['value'])
    # data_set['value'] = data_set["value"]

    #filter on frequencyband, every value in db
    df_filter1 = data_set[800:1600]
    df_filter2 = data_set[1600:3200]
    df_filter3 = data_set[3200:6400]

    #appling pa^2/hz to db
    df_filter1["value"] = pa_to_db1(df_filter1["value"])
    df_filter2["value"] = pa_to_db2(df_filter2["value"])
    df_filter3["value"] = pa_to_db3(df_filter3["value"])

    # #converting from db to ospl
    total1 = compute_ospl(list(df_filter1['value']))
    total2 = compute_ospl(list(df_filter2['value']))
    total3 = compute_ospl(list(df_filter3['value']))

# #used to plot the frequency range
#     return df_filter1
#
# check_filter = fourier_filter("mic_data/test_fwh.54.txt")
# plt.plot(check_filter)
# plt.show()


    return total1, total2, total3

#defining the list for spl data
r_spl1 = []
r_spl2 = []
r_spl3 = []


#iterating over data files in a folder called mic_data
# directory = "mic_data2" #fwh permeable
directory = "mic_data" #fwh solid
for i in range(0, 72):
    # file = "FWH_test.fwh_sim" + str(i) + ".txt"     #fwh permeable
    file = "test_fwh." + str(i) + ".txt"            #fwh solid
    filename = os.path.join(directory, file)
    print(filename)
    r1, r2, r3 = fourier_filter(filename)
    r_spl1.append(r1)
    r_spl2.append(r2)
    r_spl3.append(r3)


#defining angles for plotting
theta = np.linspace(0, 360, len(r_spl1)).tolist()

#connecting first to last point
r_spl1.append(r_spl1[0])
r_spl2.append(r_spl2[0])
r_spl3.append(r_spl3[0])
theta.append(theta[0])

#plotting
fig = go.Figure()
fig.add_trace(go.Scatterpolar(
        r = r_spl1,
        theta = theta,
        mode = 'lines',
        name = '800-1600',
        line_color = 'blue',
    ))

fig.add_trace(go.Scatterpolar(
        r = r_spl2,
        theta = theta,
        mode = 'lines',
        name = '1600-3200',
        line_color = 'green',
    ))
fig.add_trace(go.Scatterpolar(
        r = r_spl3,
        theta = theta,
        mode = 'lines',
        name = '3200-6400',
        line_color = 'red',
    ))
fig.update_polars(radialaxis=dict(visible=True,range=[25, 65]))
fig.update_layout(
    title = 'Noise Directivity at 1m around the TE in OSPL [dB]',
    showlegend = True
)

fig.show()

