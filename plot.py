# This Python file uses the following encoding: utf-8
from operation_model import Microplastic
from operation_model import Kinetics
from operation_model import Isotherm
from data_treatment import Plastic, KineticsData

# MPs/NPs
path = ["./carboxyl/1.csv", "./carboxyl/2.csv", "./carboxyl/5.csv", "./carboxyl/10.csv",
                "./carboxyl/20.csv", "./carboxyl/40.csv", "./carboxyl/50.csv"]
c_num = 8
c = [1, 2, 5, 10, 20, 40, 50]
c_label = ["1mg/L", "2mg/L", "5mg/L", "10mg/L", "20mg/L", "40mg/L", "50mg/L"]
type = "carboxyl"
excitation_wave_length = 518
x_limit = [0, 52]
y_limit = [0, 40]
s1 = r"./Fluorescence intensity - wavelength normal distribution curve.png"
s2 = r"./Fluorescence intensity - concentration fitting straight line diagram.png"
"""
------------------------------------------------------------------------------------------------------------------------
"""
data_a = Plastic.plastic_data(path)
# 8 is 1 more than the number of points (bars) plotted, the fluorescence intensity-wavelength normal distribution curve.
Microplastic.plot_fic_sc(data_a, c_num, type, "Chinese", c_label, y_limit, save=s1)
# 8 is 1 more than the number of points (bars) plotted, fluorescence intensity-concentration fitting straight line graph, return slope A and intercept B.
a1, b1 = Microplastic.plot_fic_lc(data_a, excitation_wave_length, c_num, c, type, "Chinese", x_limit, y_limit, s2)

"""
************************************************************************************************************************
************************************************************************************************************************
"""

# Kinetics
t = [0, 10, 20, 30, 45, 60, 90, 120, 180, 240]
qt = ["./carboxyl_adsorption/0min.csv", "./carboxyl_adsorption/10min.csv",
     "./carboxyl_adsorption/20min.csv", "./carboxyl_adsorption/30min.csv",
     "./carboxyl_adsorption/45min.csv", "./carboxyl_adsorption/60min.csv",
     "./carboxyl_adsorption/90min.csv", "./carboxyl_adsorption/120min.csv",
     "./carboxyl_adsorption/180min.csv", "./carboxyl_adsorption/240min.csv"]
mass = 0.036
init_volume = 0.1
sample = 0.003
c_qe = 23
k_x_limit = [0, 250]
f_y_limit = [-0.5, 2.5]
s_y_limit = [0, 15]
s3 = r"./1.png"
s4 = r"./2.png"
"""
------------------------------------------------------------------------------------------------------------------------
"""
# y1
f = KineticsData.kinetics_pfo_y(qt, excitation_wave_length, a1, b1, mass, init_volume, sample, c_qe)
# y2
s = KineticsData.kinetics_pso_y(qt, excitation_wave_length, a1, b1, mass, init_volume, sample, t)
Kinetics.kinetics_pfo(t[1:], f, type, 1, "Chinese", k_x_limit, f_y_limit, s3)
Kinetics.kinetics_pso(t[1:], s, type, 2, "Chinese", k_x_limit, s_y_limit, s4)

"""
************************************************************************************************************************
************************************************************************************************************************
"""

# isotherm
ce = [49.85, 78.32, 101.45, 125.63, 149.88]
qe = [59.78, 95.43, 143.28, 160.22, 180.01]
K_Langmuir_predict = 0.01
qm = 60
K_Freundlich_predict = 20
n = 0.3
i_x_limit = [40, 160]
i_y_limit = [50, 200]
s5 = r"./Langmuir.png"
s6 = r"./Freundlich.png"
"""
------------------------------------------------------------------------------------------------------------------------
"""
# Langmuir
Isotherm.isotherm_l(ce, qe, qm, K_Langmuir_predict, type, "Chinese", i_x_limit, i_y_limit, s5)
# Freundlich
Isotherm.isotherm_f(ce, qe, n, K_Freundlich_predict, type, "Chinese", i_x_limit, i_y_limit, s6)
