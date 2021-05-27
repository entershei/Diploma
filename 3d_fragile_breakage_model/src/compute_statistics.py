import math
import time

import numpy as np

import parameters
from draw_plots import build_parameters_for_plot_title, draw_plots

from generate_directories_names import get_cycles_info_dir, get_experiments_dir
from utils import (
    read_experiments_cycles_info,
    generate_cycle_type,
    define_cycles_representative,
)
from pathlib import Path
from statistics import mean


# Возвращает нормированное число циклов длины m, посчитанных через аналитическую формулу, указывает число циклов по
# типам (например, AAAB-циклы) и их общее количество.
def compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha):
    p_ab = 1 - p_aa - p_bb
    beta = 1 - alpha
    eps_zero = 1e-9

    def cycles_depends_on_cnt_a(l):
        if p_ab < eps_zero or alpha < eps_zero or beta < eps_zero:
            return 0.0

        r = m - l

        return (
            (x * p_ab / (alpha * beta)) ** (m - 1)
            * (alpha * r) ** (l - 1)
            / math.factorial(r)
            * (beta * l) ** (r - 1)
            / math.factorial(l)
            * alpha
            * beta
            * (1 + 2 * beta * l * p_aa / (alpha * r * p_ab)) ** (l - 1)
            * (1 + 2 * alpha * r * p_bb / (beta * l * p_ab)) ** (r - 1)
            * math.exp(
                -x
                * (
                    2 * beta * l * p_aa
                    + alpha * r * p_ab
                    + beta * l * p_ab
                    + 2 * alpha * r * p_bb
                )
                / (alpha * beta)
            )
        )

    if alpha < eps_zero:
        all_a = 0.0
    else:
        all_a = (
            (2 * x * p_aa) ** (m - 1)
            * m ** (m - 2)
            / math.factorial(m)
            / alpha ** (m - 2)
            * math.exp(-x * m * (2 * p_aa + p_ab) / alpha)
        )

    if beta < eps_zero:
        all_b = 0.0
    else:
        all_b = (
            (2 * x * p_bb) ** (m - 1)
            * m ** (m - 2)
            / math.factorial(m)
            / beta ** (m - 2)
            * math.exp(-x * m * (2 * p_bb + p_ab) / beta)
        )

    cycles = {generate_cycle_type(m, m): all_a, generate_cycle_type(m, 0): all_b}

    sum_all = all_a + all_b
    for l in range(1, m):
        cur_cycles = cycles_depends_on_cnt_a(l)
        cycles[generate_cycle_type(m, l)] = cur_cycles
        sum_all += cur_cycles

    return {"types": cycles, "all": sum_all}


# Normalized min distance
def compute_analytically_d_n(x, p_aa, p_bb, alpha, max_m):
    dist = 1.0

    for m in range(1, max_m):
        dist -= compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]

    return dist


# For non trivial cycles normalized sum of cycles lens
def compute_analytically_b_n(x, p_aa, p_bb, alpha):
    return 1 - compute_analytical_cycles_m(1, x, p_aa, p_bb, alpha)["all"]


def compute_d_by_cycles(cycles_m):
    d = 0.0
    for c_len in cycles_m.keys():
        if int(c_len) > 1:
            d += cycles_m[c_len] * (int(c_len) - 1)
    return d


def compute_b_by_cycles(cycles_m):
    b = 0.0
    for c_len in cycles_m.keys():
        if int(c_len) > 1:
            b += cycles_m[c_len] * int(c_len)
    return b


def compute_empirical_d(graph):
    return compute_d_by_cycles(graph.cycles_m)


def compute_empirical_b(graph):
    return compute_b_by_cycles(graph.cycles_m)


def draw_analytical_c_m_b(max_m, mid_m, parameters_index):
    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[parameters_index]
    parameters_str, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    step = 1e-4
    l = step
    r = 1.5 + step
    sum_c_n = []
    c2_b = []
    sum_mid_b = []
    sum_c_b = []
    for x in np.arange(l, r, step):
        cur_sum_n = 0
        cur_sum_mid_b = 0
        cur_b_n = compute_analytically_b_n(x, p_aa, p_bb, alpha)
        for m in range(2, max_m + 1):
            c_m = compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]
            cur_sum_n += c_m
            if m <= mid_m:
                cur_sum_mid_b += c_m

        sum_mid_b.append(cur_sum_mid_b / cur_b_n)
        sum_c_n.append(cur_sum_n)
        sum_c_b.append(cur_sum_n / cur_b_n)
        c2_b.append(
            compute_analytical_cycles_m(2, x, p_aa, p_bb, alpha)["all"] / cur_b_n
        )

    draw_plots(
        np.arange(l, r, step),
        [
            {
                "plot": sum_c_b,
                "label": "(c2+... +c" + str(max_m) + ")/b",
                "color": "grey",
            },
            {"plot": c2_b, "label": "c2/b", "color": "pink"},
            {
                "plot": sum_mid_b,
                "label": "(c2+...+c" + str(mid_m) + ")/b",
                "color": "blue",
                "linestyle": "dashed",
            },
        ],
        "x",
        "(Sum c_m) / b",
        "(c2 + ... + c" + str(max_m) + ") / b computed analytically",
        "3d_fragile_breakage_model/plots/statistic/analytical_c_m_b/"
        + str(max_m)
        + "_"
        + parameters_str,
    )


def estimate_p_ab(x, p_aa, alpha, aa_cycles, aaa_cycles):
    log = math.log(aaa_cycles * alpha / (aa_cycles * 2 * x * p_aa))
    return -2 * p_aa - log * alpha / x


def estimate_alpha(x, p_aa, p_bb, aa_cycles, ab_cycles, bb_cycles):
    def check_alphas(alpha1, alpha2):
        if alpha1 <= 0 and alpha2 <= 0:
            return []
        elif 1 <= alpha1 and 1 <= alpha2:
            return []
        elif 0 < alpha1 < 1 and 0 < alpha2 < 1:
            return [alpha1, alpha2]

        if alpha1 <= 0 or 1 <= alpha1:
            alpha1, alpha2 = alpha2, alpha1
        return [alpha1]

    def estimate_alpha_aa_bb():
        if p_aa != 0 and p_bb != 0 and aa_cycles != 0 and bb_cycles != 0:
            aa_cycles_divide_bb = aa_cycles / bb_cycles
            y = math.log(aa_cycles_divide_bb * p_bb / p_aa) / x
            if abs(y) < eps_zero:
                alpha1 = (p_aa - p_bb + 1) / 2
                return check_alphas(alpha1, -1)

            d = 16 + y ** 2 + 8 * y * p_aa - 8 * y * p_bb
            if d < 0:
                return []
            if d < eps_zero:
                alpha1 = (y - 4) / (2 * y)
                return check_alphas(alpha1, -1)

            alpha1 = (math.sqrt(d) + y - 4) / (2 * y)
            alpha2 = (-1 * math.sqrt(d) + y - 4) / (2 * y)
            return check_alphas(alpha1, alpha2)
        return []

    def estimate_alpha_aa_ab():
        if p_aa != 0 and p_ab != 0 and aa_cycles != 0 and ab_cycles != 0:
            aa_cycles_divide_ab = aa_cycles / ab_cycles
            y = math.log(aa_cycles_divide_ab * p_ab / p_aa) / x
            if abs(y) < eps_zero:
                alpha1 = (p_aa - p_bb + 1) / 2
                return check_alphas(alpha1, -1)

            d = 4 + y ** 2 - 4 * y * p_bb + 4 * y * p_aa
            if d < 0:
                return []
            if d < eps_zero:
                alpha1 = (-2 + y) / (2 * y)
                return check_alphas(alpha1, -1)

            alpha1 = (math.sqrt(d) - 2 + y) / (2 * y)
            alpha2 = (-1 * math.sqrt(d) - 2 + y) / (2 * y)
            return check_alphas(alpha1, alpha2)
        return []

    def estimate_alpha_bb_ab():
        if p_ab != 0 and p_bb != 0 and ab_cycles != 0 and bb_cycles != 0:
            bb_cycles_divide_ab = bb_cycles / ab_cycles
            y = math.log(bb_cycles_divide_ab * p_ab / p_bb) / x
            if abs(y) < eps_zero:
                alpha1 = (2 * p_aa + p_ab) / 2
                return check_alphas(alpha1, -1)

            d = y ** 2 + 4 + 4 * y * p_bb - 4 * y * p_aa
            if d < 0:
                return []
            if d < eps_zero:
                alpha1 = (y + 2) / (2 * y)
                return [alpha1]

            alpha1 = (math.sqrt(d) + y + 2) / (2 * y)
            alpha2 = (-1 * math.sqrt(d) + y + 2) / (2 * y)
            return check_alphas(alpha1, alpha2)
        return []

    if aa_cycles == 0 and bb_cycles == 0:
        return []
    if aa_cycles == 0 and ab_cycles == 0:
        return [0]
    if bb_cycles == 0 and ab_cycles == 0:
        return [1]

    eps_zero = 1e-8
    p_ab = max(1 - p_aa - p_bb, 0.0)
    alpha_aa_bb = estimate_alpha_aa_bb()
    alpha_aa_ab = estimate_alpha_aa_ab()
    alpha_bb_ab = estimate_alpha_bb_ab()

    res_alpha = []
    dif_p_aa_p_bb = abs(p_aa - p_bb)
    dif_p_aa_p_ab = abs(p_aa - p_ab)
    dif_p_bb_p_ab = abs(p_bb - p_ab)

    if alpha_aa_bb:
        res_alpha = alpha_aa_bb

    if alpha_bb_ab and (not res_alpha or dif_p_bb_p_ab < dif_p_aa_p_bb):
        res_alpha = alpha_bb_ab

    if alpha_aa_ab and (
        not res_alpha
        or (dif_p_aa_p_ab < dif_p_aa_p_bb and dif_p_aa_p_ab < dif_p_bb_p_ab)
    ):
        res_alpha = alpha_aa_ab

    return res_alpha


def estimate_alphas_for_graphs(start_ind, end_ind, to_represent):
    for parameter in parameters.PROBABILITIES_WITH_ALPHA[start_ind:end_ind]:
        parameters_str, p_aa, p_bb, alpha = (
            parameter["parameters_str"],
            parameter["p_aa"],
            parameter["p_bb"],
            parameter["alpha"],
        )
        print(parameters_str)
        file = (
            get_cycles_info_dir(
                parameter["number_of_experiments"],
                parameter["experiments_in_one_bunch"],
            )
            + parameters_str
            + ".csv"
        )

        Path("3d_fragile_breakage_model/plots/statistic/computed_alpha/").mkdir(
            parents=True, exist_ok=True
        )

        graphs = read_experiments_cycles_info(
            file, 4, 4, False, is_cycles_ordered=False
        )[0][:1501]

        alphas = []
        n = parameters.NUMBER_OF_FRAGILE_EDGES
        for k, graph in enumerate(graphs):
            if k == 0:
                continue

            cur_alpha = estimate_alpha(
                k / n,
                p_aa,
                p_bb,
                graph.cycles_with_edges_order[to_represent["AA"]],
                graph.cycles_with_edges_order[to_represent["AB"]],
                graph.define_cycles_representative[to_represent["BB"]],
            )
            if len(cur_alpha) > 1:
                print("two alphas")
            if cur_alpha:
                alphas.append(cur_alpha[0])

        steps = len(graphs)
        parameters_for_plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
        save_path = (
            "3d_fragile_breakage_model/plots/statistic/computed_alpha/" + parameters_str
        )

        draw_plots(
            range(1, steps),
            [
                {
                    "plot": alphas,
                    "label": "Estimated alpha",
                    "color": "blue",
                },
                {
                    "plot": [alpha] * (steps - 1),
                    "label": "Empirical alpha",
                    "color": "red",
                    "linestyle": "dashed",
                },
            ],
            "steps",
            "alpha",
            "Computed alpha by using p_aa, p_bb and number of AA, AB and BB cycles\n"
            + parameters_for_plot_title,
            save_path,
        )


def cycles_distribution_in_one_length(cycle_types, parameter_index, colors):
    start_time = time.time()
    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    to_represent, _ = define_cycles_representative(len(cycle_types[0]) + 1)
    cycle_types = list(map(lambda c_type: to_represent[c_type], cycle_types))
    experiments = read_experiments_cycles_info(
        get_experiments_dir(parameter["experiments_in_one_bunch"])
        + parameter["parameters_str"]
        + ".csv",
        len(cycle_types[0]) + 1,
        len(cycle_types[0]) + 1,
        is_int=False,
        is_cycles_ordered=True,
    )
    print("finish read", (time.time() - start_time) / 60, "m")

    mean_cycle_types = {}
    for cycle_type in cycle_types:
        mean_cycle_types[cycle_type] = []

    steps = len(experiments[0])
    for k in range(steps):
        for cycle_type in cycle_types:
            mean_cycle_types[cycle_type].append(
                mean(
                    list(
                        map(
                            lambda experiment: experiment[k].cycles_with_edges_order[
                                cycle_type
                            ],
                            experiments,
                        )
                    )
                )
            )
    print("finish counting mean", (time.time() - start_time) / 60, "m")

    cycle_types_x = {}
    xs = []
    for cycle_type in cycle_types:
        cycle_types_x[cycle_type] = []

    m = str(steps)
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    # first_divide_second = []

    for k in range(steps):
        xs.append(k / n)
        # if mean_cycle_types[cycle_types[1]][k] == 0.0:
        #     first_divide_second.append(0.0)
        # else:
        #     first_divide_second.append(mean_cycle_types[cycle_types[0]][k] / mean_cycle_types[cycle_types[1]][k])

        s = 0.0
        for cycle_type in cycle_types:
            s += mean_cycle_types[cycle_type][k]

        for cycle_type in cycle_types:
            if s == 0.0:
                cycle_types_x[cycle_type].append(0.0)
            else:
                cycle_types_x[cycle_type].append(mean_cycle_types[cycle_type][k] / s)

    # plots = [{"plot": first_divide_second, "label": cycle_types[0] + "/" + cycle_types[1], "color": colors[0]}]

    plots = []
    for i, cycle_type in enumerate(cycle_types):
        plots.append(
            {
                "plot": cycle_types_x[cycle_type],
                "label": cycle_type + "/ sum_types",
                "color": colors[i],
            }
        )

    cycle_types_str = ""
    for i, cycle_type in enumerate(cycle_types):
        cycle_types_str += cycle_type
        if i != len(cycle_types) - 1:
            cycle_types_str += ","

    # cycle_types_str = cycle_types[0] + "_" + cycle_types[1]

    save_as = (
        "3d_fragile_breakage_model/plots/statistic/type_distribution/"
        + cycle_types_str
        + "_"
        + parameter["parameters_str"]
    )

    draw_plots(
        xs,
        plots,
        "x",
        "cycles",
        "Normalized cycles for length "
        + m
        + "\n"
        + build_parameters_for_plot_title(
            parameter["p_aa"], parameter["p_bb"], parameter["alpha"]
        ),
        save_as,
    )


def main():
    colors = ["red", "blue"]
    cycles_distribution_in_one_length(["AABB", "ABAB"], 15, colors)
    cycles_distribution_in_one_length(["BBAAA", "BABAA"], 15, colors)
    cycles_distribution_in_one_length(["BBBAA", "BBABA"], 15, colors)

    cycles_distribution_in_one_length(["AABB", "ABAB"], 14, colors)
    cycles_distribution_in_one_length(["BBAAA", "BABAA"], 14, colors)
    cycles_distribution_in_one_length(["BBBAA", "BBABA"], 14, colors)

    # draw_analytical_c_m_b(max_m=10, mid_m=5, parameters_index=4)
    # draw_analytical_c_m_b(max_m=15, mid_m=10, parameters_index=4)
    # draw_analytical_c_m_b(max_m=20, mid_m=10, parameters_index=4)
    # draw_analytical_c_m_b(max_m=40, mid_m=20, parameters_index=4)
    # estimate_alphas_for_graphs(0, 5)
    # estimate_alphas_for_graphs(-4, -1)


if __name__ == "__main__":
    main()
