# This Python file uses the following encoding: utf-8
import math
import pandas as pd
import numpy as np
from copy import deepcopy as dcp
from math import log

water = pd.read_csv("water.csv").iloc[20::, 1]
wave_arr = np.array(pd.read_csv("water.csv").iloc[20::, 0], dtype=np.double)


class Plastic(object):
    # Unify the class attribute wavelength format to DataFrame
    fic_table = pd.DataFrame(wave_arr)

    @classmethod
    def plastic_data(cls, kind=None, save=None, path=None):
        """
        :description: 荧光强度-波长正态分布数据
        :param kind: MPs/NPs的改性类型
        :param save: excel数据保存路径
        :param path: 将从荧光仪取得的数据文件(.csv)路径进行按顺序存储进列表path传入
        :return:
        """
        # 防止对类属性进行叠加，对其进行深拷进行使用
        new_fic_table = dcp(cls.fic_table)

        # fic用来存储后续差分矩阵得到的新的每一列即每个浓度的新数据
        fic = []

        if path is None:
            plastic = []
        else:
            plastic = path

        # 先将水的荧光强度与每种浓度的荧光强度合并为一个表
        # 对表进行类型转换为矩阵，且数据类型是double
        # 利用差分矩阵直接后者减前者得到新的矩阵数据后再转回DataFrame
        # 将新的DataFrame的数据存入列表fic
        for temp in plastic:
            data = pd.concat([water, pd.read_csv(temp).iloc[20::, 1]],
                             axis=1, join="outer", ignore_index=True)
            arr = np.array(data, dtype=np.double)
            new_arr = np.diff(arr)
            table = pd.DataFrame(new_arr)
            fic.append(table)

        # 波长列与每种浓度的数据列整合为一个大table
        for i in range(len(fic)):
            new_fic_table = pd.concat([new_fic_table, fic[i]],
                                      axis=1,  join="outer", ignore_index=True)
        if kind:
            if save:
                save_path = "".join(["D:\\ZGJ\\SR_DATA\\", save, "\\%s\\%s.xlsx" % (kind, kind)])
            # 将new_fic_table以excel的形式保存到下列路径中
                new_fic_table.to_excel(save_path)
        # 返回这个table
        return new_fic_table

"""
*****************************************************************************************************************************************************
以下为某一浓度的微纳塑料溶液开始吸附实验，故path的取值与上述不同，是时间梯度下的荧光强度数据
"""
class Adsorption(object):
    @classmethod
    def concentration(cls, path, excitation_wave_length, a, b):
        """
        :description: 某一浓度的微纳塑料溶液在时间梯度下的吸附实验的浓度变化数据运算
        :param path: 将从荧光仪取得的数据文件(.csv)路径进行按顺序存储进列表path传入
        :param excitation_wave_length: 微纳塑料的激发波长
        :param a: 由operation_model中绘制出的荧光强度-浓度拟合直线得到的斜率a
        :param b: 由operation_model中绘制出的荧光强度-浓度拟合直线得到的截距b
        :return:
        """

        # 仿照Plastic对path传进来的数据作差分矩阵运算
        # 此时整个data表的行索引依旧是波长，而列索引变成了时间梯度
        data = Plastic.plastic_data(kind=None, save=None, path=path)
        wave_length = list(data.iloc[::, 0])
        y = []

        # 以激发波长为准，选取激发波长处的荧光强度
        # 将时间梯度中的每种荧光强度存储进y
        for val in wave_length:
            if val == excitation_wave_length:
                val_index = wave_length.index(val)
                for i in range(len(path)):
                    y.append(data.iloc[val_index, i + 1])
                break
        cnc = []
        # 利用拟合直线y=ax+b算出时间梯度中微纳塑料的每种浓度，并存入cnc列表
        for round in y:
            cnc.append((round - b) / a)
        return cnc


    @classmethod
    def adsorption_quantity(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, kind=None, concentration=None, save=None):
        """
        :description: 从0时刻开始计算时间梯度下吸附量的变化数据
        :param path: 将从荧光仪取得的数据文件(.csv)路径进行按顺序存储进列表path传入
        :param excitation_wave_length: 微纳塑料的激发波长
        :param a: 由operation_model中绘制出的荧光强度-浓度拟合直线得到的斜率a
        :param b: 由operation_model中绘制出的荧光强度-浓度拟合直线得到的截距b
        :param mass: 吸附剂气凝胶的质量
        :param init_volume: 微纳塑料溶液的初始体积
        :param sample: 时间梯度取出的小量微纳塑料样品用于测定每时刻的荧光强度
        :return:
        """

        # 得到锥形瓶中每时刻的微纳塑料溶液浓度
        cnc = cls.concentration(path, excitation_wave_length, a, b)

        # 0时刻锥形瓶中微纳塑料的总质量
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

        if kind:
            if save:
                save_path = "".join(["D:\\ZGJ\\SR_DATA\\", save, "\\%s\\%s_quantity.xlsx" % (kind, concentration)])
            # 将adsorption_quantity以excel的形式保存到下列路径中
                quantity_table = pd.DataFrame(quantity)
                quantity_table.to_excel(save_path)
        return quantity

"""
*****************************************************************************************************************************************************
以下为吸附动力学，是在某一浓度下吸附量与时间梯度的关系，故path依旧是时间梯度下的数据
"""
class KineticsData(object):
    @classmethod
    def kinetics_pfo_y(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, qe, kind=None, concentration=None, save=None):
        # 利用时间梯度下的吸附量计算动力学一阶方程拟合直线需要的y值
        qt = Adsorption.adsorption_quantity(path, excitation_wave_length, a, b, mass, init_volume, sample, kind, concentration, save)
        y = []
        e = math.e
        for i in range(len(qt)):
            if i == 0:
                continue
            y.append(log(qe - qt[i], e))
        if kind:
            if save:
                save_path = "".join(["D:\\ZGJ\\SR_DATA\\", save, "\\%s\\%s_pfo_y.xlsx" % (kind, concentration)])
            # 将kinetics_pfo_y以excel的形式保存到下列路径中
                y_table = pd.DataFrame(y)
                y_table.to_excel(save_path)
        return y

    @classmethod
    def kinetics_pso_y(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, t, kind=None, concentration=None, save=None):
        # 利用时间梯度下的吸附量计算动力学二阶方程拟合直线需要的y值
        qt = Adsorption.adsorption_quantity(path, excitation_wave_length, a, b, mass, init_volume, sample, kind, concentration, save)
        y = []
        e = math.e
        for i in range(len(qt)):
            if i == 0:
                continue
            y.append(t[i] / qt[i])
        if kind:
            if save:
                save_path = "".join(["D:\\ZGJ\\SR_DATA\\", save, "\\%s\\%s_pso_y.xlsx" % (kind, concentration)])
            # 将kinetics_pso_y以excel的形式保存到下列路径中
                y_table = pd.DataFrame(y)
                y_table.to_excel(save_path)
        return y

"""
*****************************************************************************************************************************************************
以下为吸附等温线，是在浓度梯度下饱和吸附量的变化关系，故path应取每种浓度下
上述时间梯度中得到的饱和吸附量
"""
class IsothermData(object):
    qe = []

    @classmethod
    def isotherm_l_y(cls, path, excitation_wave_length, a, b, mass, init_volume, sample, kind=None, save=None):
        for i in range(len(path)):
            temp = Adsorption.adsorption_quantity(path[i], excitation_wave_length, a, b, mass, init_volume, sample,
                                           kind=None, concentration=None, save=None)
            temp.sort()
            cls.qe.append(temp[len(temp) - 1] + 1)
        if kind:
            if save:
                save_path = "".join(["D:\\ZGJ\\SR_DATA\\", save, "\\%s\\%s_isotherm_y.xlsx" % (kind, kind)])
            # 将isotherm_l_y以excel的形式保存到下列路径中
                y_table = pd.DataFrame(cls.qe)
                y_table.to_excel(save_path)
        return cls.qe



if __name__ == "__main__":
    pass

