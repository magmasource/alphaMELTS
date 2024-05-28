from meltsengine import MELTSengine
from meltsdynamic import MELTSdynamic

import matplotlib.pyplot as plt
import numpy as np
#from scipy import stats
import math
from tasplot import *

class MELTSplotter(object):

    def demo(self, list):
        self.plotTAS(list)
        """
        self.harkerPlot(list, 'SiO2', ['CaO'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'K2O'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2', 'Na2O'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2', 'Na2O', 'FeO'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2', 'Na2O', 'FeO', 'MnO'], True)
        self.harkerPlot(list, 'SiO2', ['CaO', 'MgO', 'K2O', 'TiO2', 'Na2O', 'FeO', 'MnO', 'Al2O3'], True)
        """


    def plotTAS(self, list):
        x = list.getListProperty('dispComposition', ['liquid1'], 'sio2')
        k2o = list.getListProperty('dispComposition', ['liquid1'], 'k2o')
        na2o = list.getListProperty('dispComposition', ['liquid1'], 'na2o')
        y = [k + n for k, n in zip(k2o, na2o)]

        x,y = zip(* filter( lambda z: all(v is not None for v in z), zip(x,y) ))

        fig, axes = plt.subplots()
        add_LeMaitre_fields(axes)
        plt.scatter(x, y, marker='*', color='black')
        plt.show()

    def harkerPlot(self, list, xAxis, yAxes, plotBestFit):
        x = list.getListProperty('dispComposition', ['liquid1'], xAxis.lower())
        ys = [list.getListProperty('dispComposition', ['liquid1'], yAxis.lower()) for yAxis in yAxes]

        n = len(yAxes)

        if n == 0:
            pass

        else:
            fig, axes = plt.subplots(nrows = math.ceil(n / 2), ncols = min(n, 2), figsize=(3.6 * min(n, 2), 2.4 + 1.2 * math.ceil(n / 2)))
            plt.tight_layout(pad=3, w_pad=3, h_pad=2)

            if n != 1 and n % 2 == 1:
                axes[n // 2, 1].remove()

            for i in range(n):
                yAxis = yAxes[i]
                if n == 1:
                    axis = axes
                elif n == 2:
                    axis = axes[i]
                else:
                    axis = axes[i // 2, i % 2]
                axis.set_xlabel("{}%".format(xAxis))
                axis.set_ylabel("{}%".format(yAxis))
                axis.scatter(x, ys[i], marker='.', color='black')

                #if not plotBestFit:
                #    continue

                #slope, intercept, r_value, p_value, std_err = stats.linregress(x, ys[i])
                #line = [slope * p +intercept for p in x]
                # line = slope * x + intercept
                #axis.plot(x, line, linestyle='--')


        plt.show()