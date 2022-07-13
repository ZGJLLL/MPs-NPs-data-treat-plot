# This Python file uses the following encoding: utf-8
import numpy as np
from operation_model import Microplastic
from operation_model import Kinetics
from operation_model import Isotherm
from data_treatment import Plastic, KineticsData, IsothermData
from data_treatment import Adsorption

grade = "M_1/2st"
type = "amino"
# MPs/NPs
c = [5, 10, 25, 50, 80, 85, 90, 100, 120, 150]
c_num = 11
excitation_wave_length = 526
x_limit = [0, 152]
y_limit = [0, 140]
# data_save_path = "M_1\\2st\\PSNPs\\20220702"
data_save_path = None
# s1 = "".join([type, "_png/", grade,"/荧光强度-波长正态分布曲线图.png"])
# s2 = "".join([type, "_png/", grade,"/荧光强度-浓度拟合直线图.png"])
s1 = None
s2 = None
path = []
c_label = []
for ind, val in enumerate(c):
    path.append("".join([type, "/", str(val), ".csv"]))
for i in c:
    c_label.append("".join([str(i), "mg/L"]))
a1 = 0
b1 = 0
"""
------------------------------------------------------------------------------------------------------------------------
"""
try:
    # 导入浓度梯度下的微纳塑料母液荧光强度数据, 并完成差分矩阵运算处理
    data_a = Plastic.plastic_data(type, data_save_path, path)
    # 8要比绘制的点(条)数多1, 荧光强度-波长正态分布曲线图
    Microplastic.plot_fic_sc(data_a, c_num, type, "English", c_label, y_limit, s1)
    # 8要比绘制的点(条)数多1, 荧光强度-浓度拟合直线图，返回斜率a和截距b
    a1, b1 = Microplastic.plot_fic_lc(data_a, excitation_wave_length, c_num, c, type, "English", x_limit, y_limit, s2)
except Exception as result:
    print("请检查纳米塑料数据！%s" % result)
"""
************************************************************************************************************************
************************************************************************************************************************
"""
# 动力学(Kinetics)
t = [0, 5, 10, 15, 20, 30, 40, 50, 70, 90]
t_ipd = np.around(np.sqrt(t), 1)
mass = 0.02
init_volume = 0.1
sample = 0.004
concentration = "70"
c_qe = -1.9  # 每种浓度下的饱和吸附量
# 一/二阶模型的x和y坐标定义域
k_x_limit = [0, t[len(t) - 1] + 10]
f_y_limit = [-5, 5]
s_y_limit = [-6, 0]
# 粒子内扩散的x和y坐标定义域
i_x_limit = [t_ipd[0], t_ipd[len(t_ipd) - 1] + 2]
i_y_limit = [-25, 1]
# data_save_path = "M_1\\2st\\adsorption\\20220702"
data_save_path = None
# s3 = "".join([type, "_png/", grade,"/吸附动力学一阶模型.png"])
# s4 = "".join([type, "_png/", grade,"/吸附动力学二阶模型.png"])
# s5 = "".join([type, "_png/", grade,"/吸附动力学粒子内扩散模型.png"])
s3 = None
s4 = None
s5 = None
qt = []
for i in t[1:]:
    qt.append("".join([type, "_adsorption/", concentration, "-", str(i), ".csv"]))
"""
------------------------------------------------------------------------------------------------------------------------
"""

# 一阶方程的y值
f = KineticsData.kinetics_pfo_y(int(concentration), qt, excitation_wave_length, a1, b1, mass, init_volume, sample, c_qe, type, concentration, data_save_path)
# 二阶方程的y值
s = KineticsData.kinetics_pso_y(int(concentration), qt, excitation_wave_length, a1, b1, mass, init_volume, sample, t, type, concentration, data_save_path)
# 粒子内扩散的y值
ipd = Adsorption.adsorption_quantity(int(concentration), qt, excitation_wave_length, a1, b1, mass, init_volume, sample, type, concentration, data_save_path)
# 开始绘制动力学一阶拟合直线
Kinetics.kinetics_pfo(t[1:], f, type, 1, "English", k_x_limit, f_y_limit, s3)
# 开始绘制动力学二阶拟合直线
Kinetics.kinetics_pso(t[1:], s, type, 2, "English", k_x_limit, s_y_limit, s4)
# 开始绘制动力学粒子内扩散拟合直线
Kinetics.kinetics_ipd(t_ipd[1:], ipd[1:], type, "English", i_x_limit, i_y_limit, s5)
# except Exception as result:
#     print("请检查data_treatment中是否形参都对应或是否开始吸附实验")
"""
************************************************************************************************************************
************************************************************************************************************************
"""
# 等温线(isotherm)
try:
    # data_save_path = "M_1\\2st\\adsorption\\20220702"
    data_save_path = None
    path = []
    for i in range(len(c)):
        qt = []
        for j in t:
            qt.append("".join([type, "_adsorption/%d_" % c[i], str(j), ".csv"]))
        path.append(qt)
    ce = [49.85, 67.98]
    qe = IsothermData.isotherm_l_y(path, excitation_wave_length, a1, b1, mass, init_volume, sample, type, data_save_path)
    K_Langmuir_predict = 0.01
    qm = 60
    K_Freundlich_predict = 20
    n = 0.3
    i_x_limit = [40, 160]
    i_y_limit = [50, 200]
    # s5 = "".join([type, "_png/", grade,"/Langmuir模型.png"])
    # s6 = "".join([type, "_png/", grade,"/Freundlich模型.png"])
    s5 = None
    s6 = None
    # Langmuir曲线绘制
    Isotherm.isotherm_l(ce, qe, qm, K_Langmuir_predict, type, "English", i_x_limit, i_y_limit, s5)
    # Freundlich曲线绘制
    Isotherm.isotherm_f(ce, qe, n, K_Freundlich_predict, type, "English", i_x_limit, i_y_limit, s6)
except Exception as result:
    print("未知错误 %s" % result)
