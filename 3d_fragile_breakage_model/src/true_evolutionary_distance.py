import math
import time

import parameters
from generate_directories_names import get_cycles_info_dir
from compute_statistics import compute_analytically_d_n, compute_analytically_b_n
from draw_plots import draw_plots, build_parameters_for_plot_title
from utils import read_experiments_cycles_info, log_dictionaries, read_logs
import numpy as np


def compute_alpha(graph):
    return graph.a_in_non_trivial_cycles / (
        graph.a_in_non_trivial_cycles + graph.b_in_non_trivial_cycles
    )


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
    estimated_true_dist = best_analytical_x * analytical_n

    return {
        "estimated_true_dist": estimated_true_dist,
        "empirical_min_dist": empirical_d,
    }


def find_true_evolution_dist_fbm(graph):
    def compute_fbm_b_n(possible_gamma):
        return 1 - math.exp(-1 * possible_gamma)

    def d_divide_by_b(possible_gamma):
        max_m = 7

        analytical_min_d_n = 1

        for m in range(1, max_m):
            analytical_min_d_n -= (
                math.exp(-1 * possible_gamma * m)
                * (possible_gamma * m) ** (m - 1)
                / math.factorial(m)
                / m
            )

        return analytical_min_d_n / compute_fbm_b_n(possible_gamma)

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

    analytical_b_n = compute_fbm_b_n(best_analytical_x)
    analytical_n = empirical_b / analytical_b_n
    estimated_true_dist = best_analytical_x * analytical_n

    return {
        "estimated_true_dist": estimated_true_dist,
    }


def compute_true_evolutionary_distance(fixed_parameters, find_parameters, fbm):
    max_cycle_len_with_types = 6

    start_time = time.time()

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[5:6]:
        file, p_aa, p_bb, alpha = parameter

        print(file)

        graphs = read_experiments_cycles_info(
            get_cycles_info_dir() + file + ".csv",
            max_cycle_len_with_types,
            parameters.MAX_POSSIBLE_CYCLES_LEN,
            False,
        )[0]

        dist_info_fixed_parameters = []
        dist_info_found_parameters = []
        dist_info_fbm = []

        # Go through each step.
        for k, graph in enumerate(graphs):
            if k < 1:
                continue

            if k % 100 == 0:
                print(k)

            if fixed_parameters:
                dist_fixed = find_true_evolution_dist(
                    graph, p_aa, p_bb, alpha, max_m=40
                )
                dist_fixed["empirical_true_dist"] = k
                dist_info_fixed_parameters.append(dist_fixed)

            if find_parameters:
                found_alpha = compute_alpha(graph)
                found_p_aa, found_p_bb = compute_probabilities(graph)
                dist_found = find_true_evolution_dist(
                    graph, found_p_aa, found_p_bb, found_alpha, max_m=40
                )
                dist_found["empirical_true_dist"] = k
                dist_found["alpha"] = found_alpha
                dist_found["p_aa"] = found_p_aa
                dist_found["p_bb"] = found_p_bb
                dist_info_found_parameters.append(dist_found)

            if fbm:
                dist_fbm = find_true_evolution_dist_fbm(graph)
                dist_fbm["gamma"] = 2 * k / parameters.NUMBER_OF_FRAGILE_EDGES
                dist_info_fbm.append(dist_fbm)

        if fixed_parameters:
            log_dictionaries(
                dist_info_fixed_parameters,
                "3d_fragile_breakage_model/logs/true_evolution_distance_fixed/"
                + file
                + ".csv",
            )
        if find_parameters:
            log_dictionaries(
                dist_info_found_parameters,
                "3d_fragile_breakage_model/logs/true_evolution_distance_found_parameters/"
                + file
                + ".csv",
            )
        if fbm:
            log_dictionaries(
                dist_info_fbm,
                "3d_fragile_breakage_model/logs/true_evolution_distance_fbm/"
                + file
                + ".csv",
            )

        print("time: ", (time.time() - start_time) / 60, " m.")


def draw_dists(
    to_draw_main,
    to_draw_additional,
    additional_label,
    parameters_for_plot_title,
    save_as,
):
    xs = list(map(lambda v: int(v["empirical_true_dist"]), to_draw_main))
    empirical_true_distances = xs
    empirical_min_distances = list(
        map(lambda v: float(v["empirical_min_dist"]), to_draw_main)
    )
    estimated_true_distances = list(
        map(lambda v: float(v["estimated_true_dist"]), to_draw_main)
    )

    plots = [
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
    ]

    if to_draw_additional is not None:
        estimated_with_found_parameters = list(
            map(lambda v: float(v["estimated_true_dist"]), to_draw_additional)
        )
        plots.append(
            {
                "plot": estimated_with_found_parameters,
                "label": additional_label,
                "color": "green",
            }
        )

    draw_plots(
        xs=xs,
        plots=plots,
        x_label="Real number of rearrangements",
        y_label="Estimated number of rearrangements",
        title="Evolutionary distance\n" + parameters_for_plot_title,
        save_as=save_as,
    )


def draw_true_dist_for_parameters(f_name):
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA[5:6]:
        folder_name, p_aa, p_bb, alpha = cur_parameters
        plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
        print(plot_title)

        dist_info = read_logs(
            "3d_fragile_breakage_model/logs/" + f_name + "/" + folder_name + ".csv"
        )

        draw_dists(
            dist_info,
            None,
            None,
            plot_title,
            "3d_fragile_breakage_model/plots/" + f_name + "/" + folder_name,
        )


def draw_true_dist_with_additional_plot(
    additional_folder, additional_label, folder_to_save
):
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA[5:6]:
        parameters_str, p_aa, p_bb, alpha = cur_parameters
        plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
        print(plot_title)

        dist_info_fixed = read_logs(
            "3d_fragile_breakage_model/logs/true_evolution_distance_fixed/"
            + parameters_str
            + ".csv"
        )

        dist_info_found = read_logs(
            "3d_fragile_breakage_model/logs/"
            + additional_folder
            + parameters_str
            + ".csv"
        )

        draw_dists(
            dist_info_fixed,
            dist_info_found,
            additional_label,
            plot_title,
            "3d_fragile_breakage_model/plots/" + folder_to_save + parameters_str,
        )


if __name__ == "__main__":
    compute_true_evolutionary_distance(
        fixed_parameters=True, find_parameters=True, fbm=True
    )
    # draw_true_dist_for_parameters("true_evolution_distance_fixed")
    # draw_true_dist_for_parameters("true_evolution_distance_found_parameters")
    # draw_true_dist_with_additional_plot(
    #     "true_evolution_distance_found_parameters/",
    #     "Estimated true distance with found parameters",
    #     "true_evolution_distance_fixed_and_found_parameters/",
    # )
    # draw_true_dist_with_additional_plot(
    #     "true_evolution_distance_fbm/",
    #     "Estimated true distance under fbm",
    #     "true_evolution_distance_fixed_and_fbm/",
    # )
