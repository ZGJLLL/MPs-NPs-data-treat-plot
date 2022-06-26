# This Python file uses the following encoding: utf-8
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
from sklearn.metrics import r2_score
from time import localtime
from time import strftime
from scipy import optimize as op


def langmuir(Ce, qm, Ka):
    return qm * Ka * Ce / (1 + Ka * Ce)

def freundlich(Ce, n, Ka):
    return Ka * Ce ** n


class Microplastic(object):
    # Instantiate a linear fitting converter.
    line_model = LinearRegression()
    title_carboxyl_sc = {"English": "The fluorescence intensity/concentration standard curve of PS-COOH",
                   "Chinese": "PS-COOH的荧光强度-浓度标准曲线"}
    title_carboxyl_lc = {"English": "The fluorescence intensity/concentration linearity curve of PS-COOH",
                         "Chinese": "PS-COOH的荧光强度-浓度线性拟合曲线"}
    title_amino_sc = {"English": "The fluorescence intensity/concentration standard curve of PS-NH$_{2}$",
                      "Chinese": "PS-NH$_{2}$的荧光强度-浓度标准曲线"}
    title_amino_lc = {"English": "The fluorescence intensity/concentration linearity curve of PS-NH$_{2}$",
                      "Chinese": "PS-NH$_{2}$的荧光强度-浓度线性拟合曲线"}
    title_list = [title_carboxyl_sc, title_amino_sc, title_carboxyl_lc, title_amino_lc]
    x_label_sc = {"English": "Wave length(nm)", "Chinese": "波长(nm)"}
    y_label_sc = {"English": "Fluorescence intensity", "Chinese": "荧光强度(10$^{4}$)"}
    x_label_lc = {"English": "Concentration(mg/L)", "Chinese": "浓度(mg/L)"}
    y_label_lc = {"English": "Fluorescence intensity", "Chinese": "荧光强度(10$^{4}$)"}

    line_style = ['solid', (0, (5, 10)), (0, (1, 10)), (0, (3, 5, 1, 5)), (0, (3, 5, 1, 5, 1, 5)), 'dashed', 'dotted']
    line_color = []

    @classmethod
    def plot_fic_sc(cls, deal_data, line_num, kind, key, legend_labels, y_limit, save=None):
        print(strftime("start --> %Y-%m-%d %H:%M:%S", localtime()))
        fig = plt.figure(figsize=(12, 10), dpi=180)
        x = deal_data.iloc[::, 0]
        _x = [i for i in range(len(x))]
        line_list = []

        i = 0
        while i < line_num:
            if i == 0:
                i += 1
                continue
            y = deal_data.iloc[::, i]
            line = plt.plot(_x, y / 10000, ls=cls.line_style[i - 1])
            line_list.append(line)
            i += 1

        _xtick_labels = [i for i in range(int(x[0]), int(x[len(x) - 1]) + 1)]
        intervals = (x[len(x) - 1] - x[0]) / 10
        plt.xticks(_x[::int(intervals)], _xtick_labels[::int(intervals)])
        if kind == "amino":
            plt.title(cls.title_list[1][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[0][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        plt.minorticks_on()
        plt.tick_params(direction='out', left=False, width=1, length=6, labelsize=22)
        plt.tick_params(which='minor', direction='out', left=False, width=1, length=3)
        plt.xlabel(cls.x_label_sc[key], fontdict={"family": "Microsoft YaHei"}, size=22)
        plt.ylabel(cls.y_label_sc[key], fontdict={"family": "Microsoft YaHei"}, size=22)

        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['bottom'].set_position(('data', 0))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)
        print(strftime("end --> %Y-%m-%d %H:%M:%S", localtime()))
        print("-------------------------------------------------------------")

        plt.ylim(y_limit[0], y_limit[1])

        plt.legend(line_list, labels=legend_labels, loc=0, handlelength=5, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else:
            plt.savefig(save)

        plt.show()

    @classmethod
    def plot_fic_lc(cls, deal_data, excitation_wave_length, dot_num, xtick, kind, key, x_limit, y_limit, save=None):
        print(strftime("start --> %Y-%m-%d %H:%M:%S", localtime()))

        fig = plt.figure(figsize=(12, 10), dpi=180)

        x_dot = xtick
        wave_length = list(deal_data.iloc[::, 0])
        y = []
        y0 = []

        for value in wave_length:
            if value == excitation_wave_length:
                index_value = wave_length.index(value)
                i = 0
                while i < dot_num:
                    if i == 0:
                        i += 1
                        continue
                    y.append(deal_data.iloc[index_value, i])
                    y0.append(deal_data.iloc[index_value, i] / 10000)
                    i += 1
                break

        x_arr = np.array(x_dot).reshape(len(x_dot), 1)
        y_arr = np.array(y).reshape(len(y), 1)
        y0_arr = np.array(y0).reshape(len(y0), 1)

        cls.line_model.fit(x_arr, y_arr)
        y_predict = cls.line_model.predict(x_arr)
        a1 = cls.line_model.coef_[0][0]
        b = cls.line_model.intercept_[0]

        cls.line_model.fit(x_arr, y0_arr)
        y0_predict = cls.line_model.predict(x_arr)

        function_str = ""
        if b < 0:
            function_str = "y = %.4fx - %.4f" % (a1, -b)
        elif b > 0:
            function_str = "y = %.4fx + %.4f" % (a1, b)

        R_square = r2_score(y_arr, y_predict)
        r2_str = "R$^{2}$ = %.4f" % (R_square)

        plt.scatter(x_dot, y0, marker='s', s=80, facecolor='none', edgecolor='#365A87')

        y_function = plt.plot(x_arr, y0_predict, c="C1")
        y_r2 = plt.plot(0, 0, c="white")

        if kind == "amino":
            plt.title(cls.title_list[3][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[2][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        _x_ticks = [i for i in x_dot]
        plt.xticks(_x_ticks)

        plt.tick_params(direction='in', left=False, width=1, length=6, labelsize=22)

        plt.ylabel(cls.y_label_lc[key], fontdict={"family": "Microsoft YaHei", "size": 22})
        plt.xlabel(cls.x_label_lc[key], fontdict={"family": "Microsoft YaHei"}, size=22)
        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', y_limit[0]))
        ax.spines['left'].set_position(('data', x_limit[0]))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)
        print(strftime("end --> %Y-%m-%d %H:%M:%S", localtime()))
        print("-------------------------------------------------------------")

        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])
        plt.tick_params(direction='in', left=False, width=1, length=6, labelsize=22)

        plt.legend([y_function, y_r2], labels=[function_str, r2_str],
                   loc=0, handlelength=1, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else :
            plt.savefig(save)

        plt.show()
        return a1, b


class Kinetics(object):
    # TODO
    title_carboxyl_f = {"English": "Experimental data and calculated pseudo-first order curves of CNF-AG for PSNPs(-) removal",
                         "Chinese": "CNF-AG去除PSNPs(-)的动力学一阶拟合数据曲线"}
    title_carboxyl_s = {"English": "Experimental data and calculated pseudo-second order curves of CNF-AG for PSNPs(-) removal",
                         "Chinese": "CNF-AG去除PSNPs(-)的动力学二阶拟合数据曲线"}
    title_amino_f = {"English": "Experimental data and calculated pseudo-first order curves of CNF-AG for PSNPs(+) removal",
                         "Chinese": "CNF-AG去除PSNPs(+)的动力学一阶拟合数据曲线"}
    title_amino_s = {"English": "Experimental data and calculated pseudo-second order curves of CNF-AG for PSNPs(+) removal",
                         "Chinese": "CNF-AG去除PSNPs(+)的动力学二阶拟合数据曲线"}
    title_list = [title_carboxyl_f, title_amino_f, title_carboxyl_s, title_amino_s]
    x_label_fs = ["t(min)"]
    y_label_fs = ["ln(q$_{e}$ - q$_{t}$)", "t/q$_{t}$"]

    line_model = LinearRegression()

    @classmethod
    def kinetics_pfo(cls, t, qt, kind, order, key, x_limit, y_limit, save=None):
        fig = plt.figure(figsize=(12, 10), dpi=180)

        t_arr = np.array(t).reshape(len(t), 1)
        qt_arr = np.array(qt).reshape(len(qt), 1)

        cls.line_model.fit(t_arr, qt_arr)
        y_predict = cls.line_model.predict(t_arr)
        a1 = cls.line_model.coef_[0][0]
        b = cls.line_model.intercept_[0]
        function_str = ""
        if b < 0:
            function_str = "y = %.4fx - %.4f" % (a1, -b)
        elif b > 0:
            function_str = "y = %.4fx + %.4f" % (a1, b)

        R_square = r2_score(qt_arr, y_predict)
        r2_str = "R$^{2}$ = %.4f" % (R_square)

        line1 = plt.scatter(t, qt, s=80, marker='s', facecolor='none', edgecolor='#365A87')

        line2 = plt.plot(t_arr, y_predict, c="C1")
        y_r2 = plt.plot(0, 0, c="white")


        if kind == "amino":
            plt.title(cls.title_list[1][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[0][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        plt.tick_params(direction='in', width=1, length=6, labelsize=22)
        _x_ticks = [i for i in t]
        plt.xticks(_x_ticks)
        plt.xlabel(cls.x_label_fs[0], fontdict={"family": "Microsoft YaHei", "size": 22})
        plt.ylabel(cls.y_label_fs[order - 1], fontdict={"family": "Microsoft YaHei", "size": 22})

        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', y_limit[0]))
        ax.spines['left'].set_position(('data', x_limit[0]))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])

        plt.legend([line2, y_r2, line1], labels=[function_str, r2_str, kind],
                   loc=0, handlelength=1, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else:
            plt.savefig(save)

        plt.show()

    @classmethod
    def kinetics_pso(cls, t, y, kind, order, key, x_limit, y_limit, save=None):
        fig = plt.figure(figsize=(12, 10), dpi=180)

        t_arr = np.array(t).reshape(len(t), 1)
        qt_arr = np.array(y).reshape(len(y), 1)


        cls.line_model.fit(t_arr, qt_arr)
        y_predict = cls.line_model.predict(t_arr)
        a1 = cls.line_model.coef_[0][0]
        b = cls.line_model.intercept_[0]
        function_str = ""
        if b < 0:
            function_str = "y = %.4fx - %.4f" % (a1, -b)
        elif b > 0:
            function_str = "y = %.4fx + %.4f" % (a1, b)

        R_square = r2_score(qt_arr, y_predict)
        r2_str = "R$^{2}$ = %.4f" % (R_square)

        line1 = plt.scatter(t, y, s=80, marker='s', facecolor='none', edgecolor='#365A87')

        line2 = plt.plot(t_arr, y_predict, c="C1")
        y_r2 = plt.plot(0, 0, c="white")

        if kind == "amino":
            plt.title(cls.title_list[3][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[2][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        plt.tick_params(direction='in', width=1, length=6, labelsize=22)
        _x_ticks = [i for i in t]
        plt.xticks(_x_ticks)
        plt.xlabel(cls.x_label_fs[0], fontdict={"family": "Microsoft YaHei", "size": 22})
        plt.ylabel(cls.y_label_fs[order - 1], fontdict={"family": "Microsoft YaHei", "size": 22})

        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', y_limit[0]))
        ax.spines['left'].set_position(('data', x_limit[0]))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])

        plt.legend([line2, y_r2, line1], labels=[function_str, r2_str, kind],
                   loc=0, handlelength=1, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else:
            plt.savefig(save)

        plt.show()

    @classmethod
    def kinetics_ipd(cls):
        pass

    @classmethod
    def kinetics_e(cls):
        pass




class Isotherm(object):
    # TODO
    title_carboxyl_l = {"English": "Experimental data and calculated Langmuir curves of CNF-AG for PSNPs(-) removal",
        "Chinese": "CNF-AG去除PSNPs(-)的Langmuir拟合数据曲线"}
    title_carboxyl_f = {"English": "Experimental data and calculated Freundlich curves of CNF-AG for PSNPs(-) removal",
        "Chinese": "CNF-AG去除PSNPs(-)的Freundlich拟合数据曲线"}
    title_amino_l = {"English": "Experimental data and calculated Langmuir curves of CNF-AG for PSNPs(+) removal",
        "Chinese": "CNF-AG去除PSNPs(+)的Langmuir拟合数据曲线"}
    title_amino_f = {"English": "Experimental data and calculated Freundlich curves of CNF-AG for PSNPs(+) removal",
        "Chinese": "CNF-AG去除PSNPs(+)的Freundlich拟合数据曲线"}
    title_list = [title_carboxyl_l, title_amino_l, title_carboxyl_f, title_amino_f]
    x_label_lf = ["C$_{e}$(mg/L)"]
    y_label_lf = ["q$_{e}$(mg/g)"]

    @classmethod
    def isotherm_l(cls, c, sq, qm_predict, Ka_predict, kind, key, x_limit, y_limit, save=None):
        fig = plt.figure(figsize=(12, 10), dpi=180)
        ce = np.array(c)
        s_quantity = np.array(sq)
        popt, pv = op.curve_fit(langmuir, ce, s_quantity, p0=[qm_predict, Ka_predict])
        function_str = "q$_{m}$ = %.4f, K$_{L}$ = %.4f" % (popt[0], popt[1])
        y_predict = langmuir(ce, popt[0], popt[1])
        ce_new = np.arange(c[0] - 0.9, c[len(c) - 1] + 0.9, 0.01)
        y_predict_plot = langmuir(ce_new, popt[0], popt[1])
        R_square = r2_score(s_quantity, y_predict)
        print(R_square)
        line1 = plt.scatter(c, sq, s=80, marker='s', facecolor='none', edgecolor='#365A87')
        line2 = plt.plot(ce_new, y_predict_plot, c="C1")
        y_r2 = plt.plot(0, 0, c="white")

        r2_score2_str = "R$^{2}$ = %.4f" % (R_square)

        if kind == "amino":
            plt.title(cls.title_list[1][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[0][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        plt.tick_params(direction='in', width=1, length=6, labelsize=22)

        plt.xlabel(cls.x_label_lf[0], fontdict={"family": "Microsoft YaHei", "size": 22})
        plt.ylabel(cls.y_label_lf[0], fontdict={"family": "Microsoft YaHei", "size": 22})

        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', y_limit[0]))
        ax.spines['left'].set_position(('data', x_limit[0]))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])

        plt.legend([line2, y_r2, line1], labels=[function_str, r2_score2_str, kind],
                   loc=0, handlelength=1, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else:
            plt.savefig(save)

        plt.show()

    @classmethod
    def isotherm_f(cls, c, sq, n_predict, Ka_predict, kind, key, x_limit, y_limit, save=None):
        fig = plt.figure(figsize=(12, 10), dpi=180)
        ce = np.array(c)
        s_quantity = np.array(sq)
        popt, pv = op.curve_fit(freundlich, ce, s_quantity, p0=[n_predict, Ka_predict])
        function_str = "n = %.4f, K$_{F}$ = %.4f" % (popt[0], popt[1])
        y_predict = freundlich(ce, popt[0], popt[1])
        ce_new = np.arange(c[0] - 0.9, c[len(c) - 1] + 0.9, 0.01)
        y_predict_plot = freundlich(ce_new, popt[0], popt[1])
        R_square = r2_score(s_quantity, y_predict)
        print(R_square)
        line1 = plt.scatter(c, sq, s=80, marker='s', facecolor='none', edgecolor='#365A87')
        line2 = plt.plot(ce_new, y_predict_plot, c="C1")
        y_r2 = plt.plot(0, 0, c="white")

        r2_score2_str = "R$^{2}$ = %.4f" % (R_square)

        if kind == "amino":
            plt.title(cls.title_list[3][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        elif kind == "carboxyl":
            plt.title(cls.title_list[2][key], fontdict={"family": "Microsoft YaHei", "size": 22})
        else:
            print("请输入正确信息！")
            return

        plt.tick_params(direction='in', width=1, length=6, labelsize=22)

        plt.xlabel(cls.x_label_lf[0], fontdict={"family": "Microsoft YaHei", "size": 22})
        plt.ylabel(cls.y_label_lf[0], fontdict={"family": "Microsoft YaHei", "size": 22})

        ax = plt.gca()
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', y_limit[0]))
        ax.spines['left'].set_position(('data', x_limit[0]))

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

        plt.xlim(x_limit[0], x_limit[1])
        plt.ylim(y_limit[0], y_limit[1])

        plt.legend([line2, y_r2, line1], labels=[function_str, r2_score2_str, kind],
                   loc=0, handlelength=1, edgecolor='black', prop={'size': 20})

        if save is None:
            pass
        else:
            plt.savefig(save)

        plt.show()

    def isotherm_dr(self):
        pass

    def isotherm_tk(self):
        pass


if __name__ == "__main__":
    Kinetics.kinetics_pfo([0.5604620847941648, 1.0288509258702965, 1.6823067670473533, 3.022408998776567, 3.4531992605665423, 4.954126163609984, 7.6011992968463975, 8.1095398320345, 14.125469619841796],
                            [10, 20, 30, 45, 60, 90, 120, 180, 240], "carboxyl", 2, "English")

    Isotherm.isotherm_l([0.17, 0.25, 0.27, 0.65, 2.34], [13.1, 19.8, 26.5, 32.5, 32.8], 40, 5, "carboxyl", "English")