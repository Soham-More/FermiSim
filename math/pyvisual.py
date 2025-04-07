import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.widgets import Button
from functools import partial

from typing import List, Dict

import re

class PyVi:
    class PyViSection:
        def __init__(self, name, param):
            self.name = name
            self.param = param
            self.iterations = {}

        def parse_line(self, iter, vector):
            self.iterations[int(iter)] = np.fromstring(vector, dtype=np.float64, sep=',')
        
        def plot_iter(self, param, k, color, label_axes, fmt):
            plt.plot(param, self.iterations[k], fmt, color=color)
            if label_axes:
                plt.xlabel(self.param)
                plt.ylabel(self.name)

    def __load_from_file(self, file):
        with open(file) as f:
            lines = f.readlines()

            GLOBAL_MODE = 0
            PARAMETER_MODE = 1
            SECTION_MODE = 2
            SECTION_MODE_VEC = 3

            state = GLOBAL_MODE

            current_section = ''

            for line in lines:
                if m := re.match(r'\[(.*?)\]', line):
                    if m.group(1) == 'Parameter':    
                        state = PARAMETER_MODE
                    elif m.group(1) == 'Section':
                        state = SECTION_MODE
                    continue
                
                if m := re.match(r'(.*?):(.*)', line):
                    self.params[m.group(1)] = np.fromstring(m.group(2), dtype=np.float64, sep=',')
                
                if m := re.match(r'\((.*?)\)\-\>\[(.*?)\]', line):
                    self.sections[m.group(1)] = PyVi.PyViSection(m.group(1), m.group(2))
                    current_section = m.group(1)
                    state = SECTION_MODE_VEC
                
                if m := re.match(r'I\[(\d+)\]\=(.*)', line):
                    self.sections[current_section].parse_line(m.group(1), m.group(2))

    def __init__(self, filename):
        self.params : Dict[str, np.ndarray] = {}
        self.sections : Dict[str, PyVi.PyViSection] = {}
        self.__load_from_file(filename)
    
    def plot_iteration(self, section_name, iteration_no, label_axes=False, plt_format='o', color=None):
        param = self.params[self.sections[section_name].param]
        self.sections[section_name].plot_iter(param, iteration_no, color, label_axes, plt_format)

    # plot all sections that depend on parameter
    def plot_common(self, parameter : str, iter_no : int, plt_format='o'):
        i = 0
        for key, section in self.sections.items():
            if section.param != parameter:
                continue
            ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[i]]
            plt.plot(self.params[parameter], section.iterations[iter_no], color=ccolor, label=key)
            i = (i + 1) % len(mcolors.TABLEAU_COLORS)
        plt.xlabel(parameter)
        plt.legend()

    # plot these sections, that depend on the same parameter
    def plot_sections(self, parameter : str, sections : List[str], iter_no):
        i = 0
        for key in sections:
            section = self.sections[key]
            if section.param != parameter:
                continue
            ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[i]]
            plt.plot(self.params[parameter], section.iterations[iter_no], color=ccolor, label=key)
            i = (i + 1) % len(mcolors.TABLEAU_COLORS)
        plt.xlabel(parameter)
        plt.legend()

    def get_low_lim(self, data):
        minval = np.min(data)
        return minval - 0.05 * np.abs(minval)
    def get_high_lim(self, data):
        maxval = np.max(data)
        return maxval + 0.05 * np.abs(maxval)

    # displays all sections 
    def display_all_sections(self):
        self.curr_i = 0
        self.curr_iter = 0

        self.keys = list(self.sections.keys())

        fig, ax = plt.subplots() 
        fig.subplots_adjust(bottom=0.22)

        key = self.keys[0]
        section = self.sections[key]
        ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[0]]
        
        p, = plt.plot(self.params[section.param], section.iterations[self.curr_iter], color=ccolor)

        ax.title.set_text(f'{section.name} at iteration={self.curr_iter}')
        ax.set_xlabel(f'{section.param}')
        ax.set_ylabel(f'{section.name}')

        def draw_nplot(self : PyVi):
            key = self.keys[self.curr_i % len(self.sections.keys())]
            section = self.sections[key]
            ccolor = mcolors.TABLEAU_COLORS[list(mcolors.TABLEAU_COLORS.keys())[self.curr_i % len(mcolors.TABLEAU_COLORS)]]
            
            p.set_xdata(self.params[section.param])
            p.set_ydata(section.iterations[self.curr_iter])
            p.set_color(ccolor)
            
            ax.relim()
            ax.autoscale_view(True,True,True)

            ax.title.set_text(f'{section.name} at iteration={self.curr_iter}') 
            ax.set_xlabel(f'{section.param}')
            ax.set_ylabel(f'{section.name}')

            plt.draw()

        def nxt_plot(self : PyVi, event):
            self.curr_i = self.curr_i + 1
            draw_nplot(self)
        
        def prev_plot(self, event):
            self.curr_i = self.curr_i - 1
            draw_nplot(self)

        def nxt_iter(self : PyVi, event):
            self.curr_iter = (self.curr_iter + 1) % len(section.iterations)
            draw_nplot(self)
        
        def prev_iter(self, event):
            self.curr_iter = (self.curr_iter - 1) % len(section.iterations)
            draw_nplot(self)
        
        prevButtonLoc = fig.add_axes([0.7, 0.05, 0.1, 0.075])
        prevButton = Button(prevButtonLoc, label="Prev", color='pink', hovercolor='tomato')
        prevButton.on_clicked(partial(prev_plot, self))

        nextButtonLoc = fig.add_axes([0.81, 0.05, 0.1, 0.075])
        nextButton = Button(nextButtonLoc, label="Next", color='pink', hovercolor='tomato')
        nextButton.on_clicked(partial(nxt_plot, self))

        prevIterButtonLoc = fig.add_axes([0.1, 0.05, 0.2, 0.075])
        prevIterButton = Button(prevIterButtonLoc, label="Prev Iter", color='pink', hovercolor='tomato')
        prevIterButton.on_clicked(partial(prev_iter, self))

        nextIterButtonLoc = fig.add_axes([0.31, 0.05, 0.2, 0.075])
        nextIterButton = Button(nextIterButtonLoc, label="Next Iter", color='pink', hovercolor='tomato')
        nextIterButton.on_clicked(partial(nxt_iter, self))

        plt.show()

