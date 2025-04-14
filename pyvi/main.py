import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

# not neat, i know
import sys
sys.path.insert(0, 'math/')

from pyvisual import PyVi

pyvi = PyVi('data.pyvi')

def plot_band_diagram():
    pyvi.plot_sections('x', ['Ec', 'Ev', 'Fermi-Level'], 0, colors=['maroon', 'maroon', 'g--'])
    plt.ylabel('eV')

def plot_carrier_conc():
    pyvi.plot_sections('x', ['n', 'p'], 0, colors=['r', 'g'])
    plt.yscale('log')
    plt.ylabel('Carrier concentration($m^{-3}$)')
    plt.xlabel('x(m)')

def plot_setup():
    pyvi.plot_sections('x', ['doping'], 0, colors=['r'])
    plt.ylabel('Doping($m^{-3}$)')
    plt.xlabel('x(m)')

plot_band_diagram()
#plot_carrier_conc()

plt.show()

