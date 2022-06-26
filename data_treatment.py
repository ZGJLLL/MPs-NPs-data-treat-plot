# This Python file uses the following encoding: utf-8
import math
import pandas as pd
import numpy as np
from copy import deepcopy as dcp
from operation_model import Microplastic
from math import log

water = pd.read_csv("water.csv").iloc[20::, 1]
wave_arr = np.array(pd.read_csv("water.csv").iloc[20::, 0], dtype=np.double)


class Plastic(object):
    # Unify the class attribute wavelength format to DataFrame.
    fic_table = pd.DataFrame(wave_arr)

    @classmethod
    def plastic_data(cls, path=None):
        """
        :description: Fluorescence intensity-wavelength normal distribution data.
        :param path: The data file (.csv) path obtained from the fluorometer is sequentially stored into the list path.
        :return:
        """

        new_fic_table = dcp(cls.fic_table)

        fic = []

        if path is None:
            plastic = []
        else:
            plastic = path

        # First, the fluorescence intensity of water was combined with that of each concentration into a table.
        # Cast table to matrix and the data type is double.
        # Subtract the former from the latter to get the new matrix data and then go back to the DataFrame.
        # Save the new DataFrame data to the list fic.
        for temp in plastic:
            data = pd.concat([water, pd.read_csv(temp).iloc[20::, 1]], axis=1, join="outer", ignore_index=True)
            arr = np.array(data, dtype=np.double)
            new_arr = np.diff(arr)
            table = pd.DataFrame(new_arr)
            fic.append(table)

        for i in range(len(fic)):
            new_fic_table = pd.concat([new_fic_table, fic[i]], axis=1,  join="outer", ignore_index=True)

        return new_fic_table

"""
*****************************************************************************************************************************************************
The following is the adsorption experiment of MPs/NPs solution at a certain concentration, 
so the value of PATH is different from the above, and is the data under time gradient.
"""
class Adsorption(object):
    @classmethod
    def concentration(cls, path, excitation_wave_length, a, b):
        """
        :description: Calculation of concentration change data of a certain concentration of micro-nano plastic solution in adsorption experiment under time gradient.
        :param path: The data file (.csv) path obtained from the fluorometer is sequentially stored into the list path.
        :param excitation_wave_length: Excitation wavelength of MPs/NPs.
        :param a: Slope A obtained from the fluorescence intensity-concentration fitting line drawn in operation_model.
        :param b: Intercept B obtained from the fluorescence intensity-concentration fitting line plotted in operation_model.
        :return:
        """

        data = Plastic.plastic_data(path)
        wave_length = list(data.iloc[::, 0])
        y = []

        for val in wave_length:
            if val == excitation_wave_length:
                val_index = wave_length.index(val)
                for i in range(len(path)):
                    y.append(data.iloc[val_index, i + 1])
                break
        cnc = []

        for round in y:
            cnc.append((round - b) / a)
        return cnc


    @classmethod
    def adsorption_quantity(cls, path, excitation_wave_length, a, b, mass, init_volume, sample):
        """
        :description: The change data of adsorption capacity under time gradient was calculated from moment 0.
        :param path: The data file (.csv) path obtained from the fluorometer is sequentially stored into the list path.
        :param excitation_wave_length: Excitation wavelength of MPs/NPs.
        :param a: Slope A obtained from the fluorescence intensity-concentration fitting line drawn in operation_model.
        :param b: Intercept B obtained from the fluorescence intensity-concentration fitting line plotted in operation_model.
        :param mass: Quality of adsorbent.
        :param init_volume: MPs/NPs init_volume.
        :param sample:
        :return:
        """

        cnc = cls.concentration(path, excitation_wave_length, a, b)

        amount = init_volume * cnc[0]
        n = 0
        sum_cnc = []
        for sum_c in cnc:
            if sum_c == cnc[0]:
                continue
            n += sum_c
            sum_cnc.append(n)

        quantity = []
        for ind, c in enumerate(cnc):
            if c == cnc[0]:
                quantity.append((amount - init_volume * c) / mass)
                continue
            if ind < 2:
                quantity.append((amount - init_volume * c) / mass)
            else:
                quantity.append((amount - init_volume * c - sample * sum_cnc[ind - 2]) / mass)
            init_volume -= sample

        return quantity

"""
*****************************************************************************************************************************************************
The following is the adsorption kinetics, which is the relationship between adsorption capacity and time gradient at a certain concentration, 
so PATH is still the data under time gradient.
"""
class KineticsData(object):
    @classmethod
    def kinetics_pfo_y(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, qe):

        qt = Adsorption.adsorption_quantity(path, excitation_wave_length, a, b, mass, init_volume, sample)
        y = []
        e = math.e
        for i in range(len(qt)):
            if i == 0:
                continue
            y.append(log(qe - qt[i], e))
        return y

    @classmethod
    def kinetics_pso_y(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, t):

        qt = Adsorption.adsorption_quantity(path, excitation_wave_length, a, b, mass, init_volume, sample)
        y = []
        e = math.e
        for i in range(len(qt)):
            if i == 0:
                continue
            y.append(t[i] / qt[i])
            print(y)
        return y


if __name__ == "__main__":
    data_a = Plastic.plastic_data(["./amino/1.csv", "./amino/2.csv", "./amino/5.csv", "./amino/10.csv",
                                   "./amino/20.csv", "./amino/40.csv", "./amino/50.csv"])

    Microplastic.plot_fic_sc(data_a, 8, "amino", "Chinese",
                             ["1mg/L", "2mg/L", "5mg/L", "10mg/L", "20mg/L", "40mg/L", "50mg/L"])

    a1, b1 = Microplastic.plot_fic_lc(data_a, 518, 8, [1, 2, 5, 10, 20, 40, 50], "amino", "Chinese")
    f = KineticsData.kinetics_pfo_y(["./amino_adsorption/0min.csv", "./amino_adsorption/10min.csv",
                                     "./amino_adsorption/20min.csv", "./amino_adsorption/30min.csv",
                                     "./amino_adsorption/45min.csv", "./amino_adsorption/60min.csv",
                                     "./amino_adsorption/90min.csv", "./amino_adsorption/120min.csv",
                                     "./amino_adsorption/180min.csv", "./amino_adsorption/240min.csv"],
                                    518, a1, b1, 0.058, 0.1, 0.003, 60)
    print(f)