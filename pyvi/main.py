import numpy as np
import matplotlib.pyplot as plt

# not neat, i know
import sys
sys.path.insert(0, 'math/')

from pyvisual import PyVi

pyvi = PyVi('data.pyvi')

pyvi.save_all_sections('figs/')

#plt.show()

