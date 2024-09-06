import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseEvent
import pandas as pd

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from grb.parsers.xrt import XRT

class Cursor:
    """
    A cross hair cursor.
    """
    def __init__(self, ax):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        # text location in axes coordinates
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            # update the line positions
            self.horizontal_line.set_ydata([y])
            self.vertical_line.set_xdata([x])
            self.text.set_text(f'x={x:1.2f}, y={y:1.2f}')
            self.ax.figure.canvas.draw()

data = pd.read_csv('events/grb231118A/data/optical_data.csv')

# filter to only skynet data
data = data[data['Source'] == 'Skynet']

# reindex
data = data.reset_index(drop=True)

source = data['Source']
mjd = data['MJD']
mag = data['Mag']
filter = data['Filter']
mag_err = data['Uncertainty']

xrt = XRT('events/grb231118A/data/xrt_data.txt')

xrt.parse()

show_plot = input('Show XRT plot? (y/n)')
if show_plot == 'y':
    xrt.plot(wt=True)

xrt_times = np.array(xrt.windowed_timing.times + xrt.photon_counting.times)
xrt_fluxes = np.array(xrt.windowed_timing.fluxes + xrt.photon_counting.fluxes)

xrt_spectral_flux = xrt_fluxes / (2.34 * 10 ** 17) * (10 ** 20)


# Map magnitudes to spectral fluxes
conversion = {
    'B': 4.063 ,
    'V': 3.636,
    'R': 3.064,
    'I': 2.416
}

optical_spectral_flux = np.zeros(len(source))
optical_spectral_flux_err = np.zeros(len(source))
for i in range(len(source)):
    optical_spectral_flux[i] = (conversion[filter[i]]) * (10 ** (mag[i] / -2.5))
    optical_spectral_flux_err[i] = abs((conversion[filter[i]]) * (10 ** ((mag[i] + mag_err[i]) / -2.5)) - optical_spectral_flux[i])


event_time = 60266.71983

time_since_event = mjd - event_time

seconds_since_event = time_since_event * 86400

# map filters to colors
colors = {
    'u\'': 'magenta',
    'g\'': 'cyan',
    'r\'': 'orange',
    'i\'': 'black',
    'z\'': 'lightskyblue',
    'q\'': 'lightgrey',
    'B': 'blue',
    'V': 'green',
    'R': 'red',
    'I': 'black',
    'white': 'darkviolet',
    'v': 'green',
    'b': 'blue',
    'u': 'magenta',
    'w1': 'mediumorchid',
}

# map sources to markers
markers = {
    'LCOGT': 's',
    'LCO': 's',
    'Skynet': 'o',
    'MeerLICHT': 'D',
    'REM': '*',
    'VLT': 'X',
    'Swift/UVOT': 'P',
}

show_plot = input('Show optical plot? (y/n)')
if show_plot == 'y':
    plt.figure()
    for source_name, source_data in data.groupby('Source'):
        for filter_name, filter_data in source_data.groupby('Filter'):
            plt.scatter(seconds_since_event[filter_data.index], mag[filter_data.index], marker=markers[source_name], color=colors[filter_name], label=source_name + ' ' + filter_name)
            plt.errorbar(seconds_since_event[filter_data.index], mag[filter_data.index], yerr=mag_err[filter_data.index], fmt=markers[source_name], color=colors[filter_name])

    plt.xscale('log')
    plt.xlabel('Time since BAT trigger (s)')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.xlim(80, 10**6)
    plt.ylim(24, 13)
    plt.title('GRB 231118A Light Curve')
    cursor = Cursor(plt.gca())
    plt.connect('motion_notify_event', cursor.on_mouse_move)
    plt.legend(loc='lower left', fontsize='small')
    plt.show()

show_plot = input('Show combo plot? (y/n)')
if show_plot == 'y':
    plt.figure()
    for source_name, source_data in data.groupby('Source'):
        for filter_name, filter_data in source_data.groupby('Filter'):
            plt.scatter(seconds_since_event[filter_data.index], optical_spectral_flux[filter_data.index], marker=markers[source_name], color=colors[filter_name], label=source_name + ' ' + filter_name)
            plt.errorbar(seconds_since_event[filter_data.index], optical_spectral_flux[filter_data.index], yerr=optical_spectral_flux_err[filter_data.index], fmt=markers[source_name], color=colors[filter_name])
    plt.scatter(xrt_times, xrt_spectral_flux, marker='x', color='cyan', label='XRT')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time since BAT trigger (s)')
    plt.ylabel('Flux')
    plt.xlim(80, 10**6)
    plt.title('GRB 231118A Light Curve')
    cursor = Cursor(plt.gca())
    plt.connect('motion_notify_event', cursor.on_mouse_move)
    plt.legend(loc='lower left', fontsize='small')
    plt.show()