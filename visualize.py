""" module for visualizing miRNA expression values """

import matplotlib.pyplot as plt
from pylab import *
import seaborn as sns
import numpy as np
import math
from ml_tools import *

""" colors from colorbrewer2.org """
#COL = ["#e0f3db", "#a8ddb5", "#43a2ca"]
COL5 = ["#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac"]
COL = [COL5[0], COL5[1], COL5[3]]

class MirnaVisualizer:
    """ create bar chart of cleavage site expression values """

    args = {
            "col_start" : COL[1],
            "col_end"   : COL[2],
            "col_shade" :  COL[0],
            "col_5p_shade" : COL[0],
            "col_3p_shade" : COL[0],
            "hist_alpha" : .7,
            "shade_alpha" : .7,
            "figsize"   : (18, 4),
            "xtick_fontsize"    : 8
    }


    def __init__(self, name, length, **kwargs):
        """Create Plots and Graphics of a given microRNA

        Args:
            name: name of the miRNA
            length: number of positions on the x-axis

        """

        self.name = name
        self.x_range = range(length)
        for key, value in kwargs.iteritems():
            if self.args.get(key) is not None:
                self.args[key] = value

        sns.set_style("ticks", {"axes.grid" : True})
        self.fig, self.ax = plt.subplots(figsize=self.args["figsize"])
        self.ax.set_title(self.name)

    def _prepare(self):
        """stuff to do before returning the graphic"""
        self.ax.legend()

    def _twin_ax(self):
        if not hasattr(self, 'curve_ax'):
            self.curve_ax = self.ax.twinx()
            self.curve_ax.grid(False)

    def show(self):
        self._prepare()
        self.fig.show()

    def get_fig(self):
        self._prepare()
        return self.fig

    def add_xticks(self, seq):
        """set xticks provided an array with labels"""
        self.ax.set_xticks(range(len(seq)))
        self.ax.set_xticklabels(seq)
        self.ax.tick_params(axis="x", which="both",
                labelsize=self.args["xtick_fontsize"])

    def add_vline(self, pos, **kwargs):
        self.ax.axvline(pos, **kwargs)

    def add_expression(self, start_bucket, end_bucket):
        x = self.x_range

        """shorter bar in front"""
        background_start = [0] * len(x)
        background_end = [0] * len(x)
        front_start = [0] * len(x)
        front_end = [0] * len(x)
        for i in x:
            if (start_bucket[i] > end_bucket[i]):
                background_start[i] = start_bucket[i]
                front_end[i] = end_bucket[i]
            else:
                background_end[i] = end_bucket[i]
                front_start[i] = start_bucket[i]
        try:
            self.ax.bar(x, background_start, color=self.args["col_start"],
                    label="start", log=1, align="center")
            self.ax.bar(x, background_end, color=self.args["col_end"],
                    label="end", log=1, align="center")
            self.ax.bar(x, front_start, color=self.args["col_start"],
                    log=1, alpha=self.args["hist_alpha"], align="center")
            self.ax.bar(x, front_end, color=self.args["col_end"],
                    log=1, alpha=self.args["hist_alpha"], align="center")
        except ValueError:
            """ if one of the lists is empty """
            pass

    def add_mature(self, pos5p, pos3p):
        """ shade the annotated cleavage sites. Provide tuple with
        start and end for each
        """
        self.ax.axvspan(pos5p[0], pos5p[1],
                alpha=self.args["shade_alpha"], color=self.args["col_5p_shade"])
        self.ax.axvspan(pos3p[0], pos3p[1],
                alpha=self.args["shade_alpha"], color=self.args["col_3p_shade"])

    def add_log_curve(self, y_values, label, **kwargs):
        self._twin_ax()
        self.add_curve(y_values, label, **kwargs)
        self.curve_ax.set_yscale("log")

    def add_curve(self, y_values, label, **kwargs):
        self._twin_ax()
        self.curve_ax.plot(self.x_range, y_values, label=label, **kwargs)

