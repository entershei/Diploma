import math
import random

import parameters
from generate_directories_names import get_cycles_info_dir
from src.compute_statistics import compute_analytically_d_n, compute_analytically_b_n
from src.draw_plots import draw_plots
from utils import read_experiments_cycles_info, log_dictionaries, read_logs
import numpy as np


def compute_alpha(graph):
    a_edges = graph.a_in_non_trivial_cycles + graph.cycle_types["A"]
    b_edges = graph.b_in_non_trivial_cycles + graph.cycle_types["B"]
    error = random.randint(-10, 10)
    a_edges += error
    b_edges += error

    division = a_edges / b_edges
    # division = graph.a_in_non_trivial_cycles / graph.b_in_non_trivial_cycles
    alpha = division / (1 + division)
    return alpha


def compute_probabilities(graph):
    aa_cycles = graph.cycle_types["AA"]
    ab_cycles = graph.cycle_types["AB"]
    bb_cycles = graph.cycle_types["BB"]
    summ = aa_cycles + ab_cycles + bb_cycles
    return aa_cycles / summ, bb_cycles / summ


def compute_probabilities2(graph):
    aaa_cycles = graph.cycle_types["AAA"]
    bbb_cycles = graph.cycle_types["BBB"]
    cycles3 = graph.cycles_m["3"]

    p_aa = math.sqrt(aaa_cycles / cycles3 / 3)
    p_bb = math.sqrt(bbb_cycles / cycles3 / 3)

    p_ab = 1 - p_aa - p_bb

    # check = abs(p_ab ** 2 + 2 * p_ab * p_aa - graph.cycle_types["AAB"] / cycles3)

    return p_aa, p_ab, p_bb


def find_true_evolution_dist(graph, p_aa, p_bb, alpha, max_m):
    def d_divide_by_b(possible_x):
        analytical_min_d = compute_analytically_d_n(
            possible_x,
            p_aa,
            p_bb,
            alpha,
            max_m,
        )
        analytical_sum_lens = compute_analytically_b_n(possible_x, p_aa, p_bb, alpha)
        return analytical_min_d / analytical_sum_lens

    step = 5e-4
    l = step
    r = 1 + step

    empirical_d = 0
    empirical_b = 0
    for cycle_len in graph.cycles_m.keys():
        if int(cycle_len) > 1:
            empirical_d += graph.cycles_m[cycle_len] * (int(cycle_len) - 1)
            empirical_b += graph.cycles_m[cycle_len] * int(cycle_len)

    empirical_d_b = empirical_d / empirical_b

    best_analytical_x = 0
    min_error = 10000
    for x in np.arange(l, r, step):
        analytical_d_b = d_divide_by_b(x)
        cur_error = abs(analytical_d_b - empirical_d_b)
        if cur_error < min_error:
            best_analytical_x, min_error = x, cur_error

    analytical_b = compute_analytically_b_n(best_analytical_x, p_aa, p_bb, alpha)
    analytical_n = empirical_b / analytical_b
    analytical_true_dist = best_analytical_x * analytical_n

    return {
        "analytical_true_dist": analytical_true_dist,
        "real_min_dist": empirical_d,
    }


def compute():
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = 30

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[5:]:
        file, p_aa, p_bb, alpha = parameter

        print(file)

        graphs = read_experiments_cycles_info(
            get_cycles_info_dir() + file + ".csv",
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            False,
        )[0]

        dist_info = []
        # Go through each step.
        for k, graph in enumerate(graphs):
            if k < 1:
                continue

            dist = find_true_evolution_dist(graph, p_aa, p_bb, alpha, max_m=30)
            dist["real_true_dist"] = k
            dist_info.append(dist)

        log_dictionaries(
            dist_info,
            "fragile_breakage_model/logs/true_evolution_distance/" + file,
        )


def draw_dists(to_draw, parameters_for_plot_name, save_as):
    xs = list(map(lambda v: int(v["real_true_dist"]), to_draw))
    empirical_true_distances = xs
    empirical_min_distances = list(map(lambda v: float(v["real_min_dist"]), to_draw))
    estimated_true_distances = list(
        map(lambda v: float(v["analytical_true_dist"]), to_draw)
    )

    draw_plots(
        xs=xs,
        plots=[
            {
                "plot": empirical_true_distances,
                "label": "Empirical true distance",
                "color": "red",
            },
            {
                "plot": estimated_true_distances,
                "label": "Estimated true distance",
                "color": "blue",
            },
            {
                "plot": empirical_min_distances,
                "label": "Empirical min distance",
                "color": "black",
            },
        ],
        x_label="Real number of rearrangements",
        y_label="Estimated number of rearrangements",
        title="Evolutionary distance\n" + parameters_for_plot_name,
        save_as=save_as,
    )


def draw_true_dist_for_parameters():
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        folder_name, p_aa, p_bb, alpha = cur_parameters
        print(folder_name)

        dist_info = read_logs(
            "fragile_breakage_model/plots/true_evolution_distance/" + folder_name
        )

        draw_dists(
            dist_info,
            folder_name,
            "fragile_breakage_model/plots/true_evolution_distance/" + folder_name,
        )


if __name__ == "__main__":
    compute()
    draw_true_dist_for_parameters()
