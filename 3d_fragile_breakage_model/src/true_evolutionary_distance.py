import time
from statistics import median

import parameters
from generate_directories_names import get_cycles_info_dir, get_experiments_dir
from compute_statistics import (
    compute_analytically_d_n,
    compute_analytically_b_n,
    estimate_alpha,
    estimate_p_ab,
    compute_empirical_d,
    compute_empirical_b,
)
from draw_plots import draw_plots, build_parameters_for_plot_title, draw
from compute_statistics import compute_analytical_cycles_m
from utils import (
    read_experiments_cycles_info,
    log_dictionaries,
    read_logs,
    define_cycles_representative,
    generate_cycle_types_representative,
)
import matplotlib.pyplot as plt
import numpy as np
import skopt


def find_true_evolution_dist_different_errors(graph, p_aa, p_bb, alpha):
    empirical_d = compute_empirical_d(graph)
    empirical_b = compute_empirical_b(graph)
    empirical_d_b = empirical_d / empirical_b

    step = 5e-4
    l = step
    r = 1.5 + step

    min_errors = [100000] * 7
    best_analytical_xs = [0] * 7
    sum_to1 = 7
    sum_to2 = 21
    m = 20

    for x in np.arange(l, r, step):
        analytical_b_n = compute_analytically_b_n(x, p_aa, p_bb, alpha)

        d_b_error = abs(
            compute_analytically_d_n(x, p_aa, p_bb, alpha, m) / analytical_b_n
            - empirical_d_b
        )
        if d_b_error < min_errors[0]:
            best_analytical_xs[0], min_errors[0] = x, d_b_error

        c2_error = abs(
            compute_analytical_cycles_m(2, x, p_aa, p_bb, alpha)["all"] / analytical_b_n
            - graph.cycles_m["2"] / empirical_b
        )
        if c2_error < min_errors[1]:
            best_analytical_xs[1], min_errors[1] = x, c2_error

        sum_c_m_empirical1 = 0
        sum_c_m_analytical1 = 0
        separate_error1 = 0
        for i in range(2, sum_to1):
            sum_c_m_empirical1 += graph.cycles_m[str(i)]
            sum_c_m_analytical1 += compute_analytical_cycles_m(i, x, p_aa, p_bb, alpha)[
                "all"
            ]
            separate_error1 += abs(
                compute_analytical_cycles_m(i, x, p_aa, p_bb, alpha)["all"]
                / analytical_b_n
                - graph.cycles_m[str(i)] / empirical_b
            )

        sum_c_m_empirical2 = sum_c_m_empirical1
        sum_c_m_analytical2 = sum_c_m_analytical1
        separate_error2 = separate_error1

        for i in range(sum_to1, sum_to2):
            sum_c_m_empirical2 += graph.cycles_m[str(i)]
            sum_c_m_analytical2 += compute_analytical_cycles_m(i, x, p_aa, p_bb, alpha)[
                "all"
            ]
            separate_error2 += abs(
                compute_analytical_cycles_m(i, x, p_aa, p_bb, alpha)["all"]
                / analytical_b_n
                - graph.cycles_m[str(i)] / empirical_b
            )

        sum_error1 = abs(
            sum_c_m_empirical1 / empirical_b - sum_c_m_analytical1 / analytical_b_n
        )
        sum_error2 = abs(
            sum_c_m_empirical2 / empirical_b - sum_c_m_analytical2 / analytical_b_n
        )
        if sum_error1 < min_errors[2]:
            best_analytical_xs[2], min_errors[2] = (
                x,
                sum_error1,
            )
        if sum_error2 < min_errors[3]:
            best_analytical_xs[3], min_errors[3] = (
                x,
                sum_error2,
            )

        if separate_error1 < min_errors[4]:
            best_analytical_xs[4], min_errors[4] = (
                x,
                separate_error1,
            )

        if separate_error2 < min_errors[5]:
            best_analytical_xs[5], min_errors[5] = (
                x,
                separate_error2,
            )

        cur_combine_error = abs(
            compute_analytically_d_n(x, p_aa, p_bb, alpha, 20) / analytical_b_n
            - empirical_d_b
        )
        cur_combine_error += separate_error1 * 3

        if cur_combine_error < min_errors[6]:
            best_analytical_xs[6], min_errors[6] = (
                x,
                cur_combine_error,
            )

    analytical_bs = list(
        map(
            lambda best_x: compute_analytically_b_n(best_x, p_aa, p_bb, alpha),
            best_analytical_xs,
        )
    )
    analytical_ns = list(map(lambda best_b: empirical_b / best_b, analytical_bs))
    estimated_true_dists = []
    for i in range(len(analytical_ns)):
        estimated_true_dists.append(best_analytical_xs[i] * analytical_ns[i])

    res = {"empirical_min_dist": empirical_d}

    for i in range(len(estimated_true_dists)):
        res["estimated_true_dist" + str(i)] = estimated_true_dists[i]
        res["best_x" + str(i)] = best_analytical_xs[i]
        res["min_error" + str(i)] = min_errors[i]

    return res


def find_true_evolution_dist_fixed_parameters(graph, p_aa, p_bb, alpha, use_b_d_in_log):
    if use_b_d_in_log:
        empirical_d = graph.d
        empirical_b = graph.b
    else:
        empirical_d = compute_empirical_d(graph)
        empirical_b = compute_empirical_b(graph)

    step = 5e-4
    l_x = step
    r_x = 1.5 + step

    while r_x - l_x > step:
        x1 = (2 * l_x + r_x) / 3
        x2 = (l_x + 2 * r_x) / 3
        error1 = compute_error(x1, p_aa, p_bb, alpha, graph, empirical_b)
        error2 = compute_error(x2, p_aa, p_bb, alpha, graph, empirical_b)
        if error1 < error2:
            r_x = x2
        else:
            l_x = x1

    best_x = (l_x + r_x) / 2
    min_error = compute_error(best_x, p_aa, p_bb, alpha, graph, empirical_b)

    analytical_b_n = compute_analytically_b_n(best_x, p_aa, p_bb, alpha)
    analytical_n = empirical_b / analytical_b_n
    estimated_true_dist = best_x * analytical_n

    return {
        "empirical_min_dist": empirical_d,
        "estimated_true_dist": estimated_true_dist,
        "best_x": best_x,
        "min_error": min_error,
    }


def compute_error(x, p_aa, p_bb, alpha, graph, empirical_b):
    sum_to = 10
    error = 0

    analytical_b_n = compute_analytically_b_n(x, p_aa, p_bb, alpha)
    for i in range(2, sum_to):
        error += abs(
            compute_analytical_cycles_m(i, x, p_aa, p_bb, alpha)["all"] / analytical_b_n
            - graph.cycles_m[str(i)] / empirical_b
        )
    return error


def find_true_evolution_dist_with_parameters_using_skopt(graph, n_calls, to_represent):
    space = [
        skopt.space.Real(0, 1.5, name="x"),
        skopt.space.Real(0, 1, name="p_aa"),
        skopt.space.Real(0, 1 / 3, name="p_ab"),
    ]

    def evaluate_error(params):
        max_score = 1e10

        x = params["x"]
        p_aa = params["p_aa"]
        p_ab = params["p_ab"]
        p_bb = 1 - p_aa - p_ab

        if p_bb < p_ab or p_aa < p_ab:
            return max_score

        alphas = estimate_alpha(
            x,
            p_aa,
            p_bb,
            graph.cycles_with_edges_order[aa_represent],
            graph.cycles_with_edges_order[ab_represent],
            graph.cycles_with_edges_order[bb_represent],
        )
        min_error = max_score
        for alpha in alphas:
            error = compute_error(x, p_aa, p_bb, alpha, graph, empirical_b)
            if error < min_error:
                min_error = error

        return min_error

    @skopt.utils.use_named_args(space)
    def objective(**params):
        return evaluate_error(params)

    empirical_b = compute_empirical_b(graph)
    results = skopt.forest_minimize(objective, space, n_calls=n_calls)
    best_params = results.x
    best_x, best_p_aa, best_p_bb = best_params[0], best_params[1], best_params[2]
    aa_represent = to_represent["AA"]
    ab_represent = to_represent["AB"]
    bb_represent = to_represent["BB"]

    best_alphas = estimate_alpha(
        best_x,
        best_p_aa,
        best_p_bb,
        graph.cycles_with_edges_order[aa_represent],
        graph.cycles_with_edges_order[ab_represent],
        graph.cycles_with_edges_order[bb_represent],
    )
    error1 = compute_error(
        best_x, best_p_aa, best_p_bb, best_alphas[0], graph, empirical_b
    )
    best_alpha = best_alphas[0]
    if len(best_alphas) > 1:
        error2 = compute_error(
            best_x, best_p_aa, best_p_bb, best_alphas[1], graph, empirical_b
        )
        if error2 < error1:
            best_alpha = best_alphas[1]
            error1 = error2

    best_analytical_b_n = compute_analytically_b_n(
        best_x, best_p_aa, best_p_bb, best_alpha
    )
    analytical_n = empirical_b / best_analytical_b_n
    estimated_true_dist = best_x * analytical_n

    return {
        "empirical_min_dist": compute_empirical_d(graph),
        "estimated_true_dist": estimated_true_dist,
        "best_x": best_x,
        "best_p_aa": best_p_aa,
        "best_p_ab": 1 - best_p_aa - best_p_bb,
        "best_p_bb": best_p_bb,
        "best_alpha": best_alpha,
        "min_error": error1,
    }


def find_true_evolution_dist_and_find_parameters0(graph, to_represent):
    empirical_d = compute_empirical_d(graph)
    empirical_b = compute_empirical_b(graph)

    step_x = 5e-4
    l_x = step_x
    r_x = 1.5 + step_x

    min_error = 1e10

    best_analytical_x = 0.0
    best_p_aa = 0.0
    best_p_ab = 0.0
    best_p_bb = 0.0
    best_alpha = 0.0
    p_step = 0.015
    aa_represent = to_represent["AA"]
    aaa_represent = to_represent["AAA"]

    for x in np.arange(l_x, r_x, step_x):
        for p_aa in np.arange(p_step, 1, p_step):
            for alpha in np.arange(p_step, 1, p_step):
                p_ab = estimate_p_ab(
                    x,
                    p_aa,
                    alpha,
                    graph.cycles_with_edges_order[aa_represent],
                    graph.cycles_with_edges_order[aaa_represent],
                )
                p_bb = max(1 - p_aa - p_ab, 0.0)
                if p_ab < 0 or p_aa < p_ab or p_bb < p_ab:
                    continue

                error = compute_error(x, p_aa, p_bb, alpha, graph, empirical_b)
                if error < min_error:
                    best_analytical_x = x
                    best_p_aa = p_aa
                    best_p_ab = p_ab
                    best_p_bb = p_bb
                    best_alpha = alpha
                    min_error = error

    analytical_b_n = compute_analytically_b_n(
        best_analytical_x, best_p_aa, best_p_bb, best_alpha
    )
    analytical_n = empirical_b / analytical_b_n
    estimated_true_dist = best_analytical_x * analytical_n

    return {
        "empirical_min_dist": empirical_d,
        "estimated_true_dist": estimated_true_dist,
        "best_x": best_analytical_x,
        "best_p_aa": best_p_aa,
        "best_p_ab": best_p_ab,
        "best_p_bb": best_p_bb,
        "best_alpha": best_alpha,
        "min_error": min_error,
    }


def find_true_evolution_dist_and_find_parameters1(graph, to_represent, use_b_d_in_log):
    def compute_error_and_find_parameters(x):
        best_p_aa_x = 0.0
        best_p_ab_x = 0.0
        best_p_bb_x = 0.0
        best_alpha_x = 0.0
        min_error_x = 1e10

        p_step = 0.01
        aa_represent = to_represent["AA"]
        ab_represent = to_represent["AB"]
        bb_represent = to_represent["BB"]

        for p_aa in np.arange(p_step, 1, p_step):
            for p_ab in np.arange(0, min(1 / 3, (1 - p_aa) / 2), p_step):
                p_bb = max(1 - p_aa - p_ab, 0.0)
                alphas = estimate_alpha(
                    x,
                    p_aa,
                    p_bb,
                    graph.cycles_with_edges_order[aa_represent],
                    graph.cycles_with_edges_order[ab_represent],
                    graph.cycles_with_edges_order[bb_represent],
                )

                for alpha in alphas:
                    error = compute_error(x, p_aa, p_bb, alpha, graph, empirical_b)
                    if error < min_error_x:
                        best_p_aa_x = p_aa
                        best_p_ab_x = p_ab
                        best_p_bb_x = p_bb
                        best_alpha_x = alpha
                        min_error_x = error

        return min_error_x, best_p_aa_x, best_p_ab_x, best_p_bb_x, best_alpha_x

    if use_b_d_in_log:
        empirical_d = graph.d
        empirical_b = graph.b
    else:
        empirical_d = compute_empirical_d(graph)
        empirical_b = compute_empirical_b(graph)

    step_x = 5e-4
    l_x = step_x
    r_x = 1.5 + step_x

    while r_x - l_x > step_x:
        x1 = (2 * l_x + r_x) / 3
        x2 = (l_x + 2 * r_x) / 3
        error1, _, _, _, _ = compute_error_and_find_parameters(x1)
        error2, _, _, _, _ = compute_error_and_find_parameters(x2)
        if error1 < error2:
            r_x = x2
        else:
            l_x = x1

    best_x = (l_x + r_x) / 2
    (
        min_error,
        best_p_aa,
        best_p_ab,
        best_p_bb,
        best_alpha,
    ) = compute_error_and_find_parameters(best_x)

    analytical_b_n = compute_analytically_b_n(best_x, best_p_aa, best_p_bb, best_alpha)
    analytical_n = empirical_b / analytical_b_n
    estimated_true_dist = best_x * analytical_n

    return {
        "empirical_min_dist": empirical_d,
        "estimated_true_dist": estimated_true_dist,
        "best_x": best_x,
        "best_p_aa": best_p_aa,
        "best_p_ab": best_p_ab,
        "best_p_bb": best_p_bb,
        "best_alpha": best_alpha,
        "min_error": min_error,
    }


def compute_true_evolutionary_distance(parameter_index, method):
    max_cycle_len_with_types = 5
    # to_represent, _ = define_cycles_representative(max_cycle_len_with_types)
    to_represent, _ = generate_cycle_types_representative(1, max_cycle_len_with_types)

    start_time = time.time()
    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    file, p_aa, p_bb, alpha = (
        parameter["parameters_str"],
        parameter["p_aa"],
        parameter["p_bb"],
        parameter["alpha"],
    )

    print(file)

    graphs = read_experiments_cycles_info(
        get_cycles_info_dir(
            parameter["number_of_experiments"], parameter["experiments_in_one_bunch"]
        )
        + file
        + ".csv",
        max_cycle_len_with_types,
        parameters.MAX_POSSIBLE_CYCLES_LEN,
        is_int=False,
        is_cycles_ordered=False,
    )[0]

    dist_info = []

    # Go through each step.
    for k, graph in enumerate(graphs):
        if k < 1 or (k % 100) != 0:
            continue

        # cur_dist_info = find_true_evolution_dist_fixed_parameters(
        #     graph, p_aa, p_bb, alpha, use_b_d_in_log=False
        # )

        cur_dist_info = find_true_evolution_dist_and_find_parameters1(
            graph, to_represent, use_b_d_in_log=False
        )

        cur_dist_info["real_alpha"] = alpha
        cur_dist_info["real_p_aa"] = p_aa
        cur_dist_info["real_p_bb"] = p_bb

        cur_dist_info["empirical_true_dist"] = k
        dist_info.append(cur_dist_info)

        if k % 100 == 0:
            print(
                "k:",
                k,
                "best_x:",
                cur_dist_info["best_x"],
                "best_p_aa:",
                cur_dist_info["best_p_aa"],
                "best_p_ab:",
                cur_dist_info["best_p_ab"],
                "best_p_bb:",
                cur_dist_info["best_p_bb"],
                "best_alpha:",
                cur_dist_info["best_alpha"],
                "min_error:",
                cur_dist_info["min_error"],
            )
            print(
                "time:",
                (time.time() - start_time) / 60,
                " m.",
            )

    log_dictionaries(
        dist_info,
        "3d_fragile_breakage_model/logs/" + method + file + ".csv",
    )

    print("time:", (time.time() - start_time) / 60, " m.")


def compute_true_evolutionary_distance_for_box_plot(
    parameter_index, number_of_experiments, k_step, method
):
    max_cycle_len_with_types = 4
    sum_to = 10
    to_represent, _ = define_cycles_representative(max_cycle_len_with_types)
    start_time = time.time()
    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    file, p_aa, p_bb, alpha = (
        parameter["parameters_str"],
        parameter["p_aa"],
        parameter["p_bb"],
        parameter["alpha"],
    )
    if parameter["experiments_in_one_bunch"] != 1:
        print("experiments_in_one_bunch should be 1")
        return

    print(file)

    experiments = read_experiments_cycles_info(
        get_experiments_dir(parameter["experiments_in_one_bunch"]) + file + ".csv",
        max_cycle_len_with_types,
        sum_to,
        is_int=False,
        is_cycles_ordered=True,
    )[:number_of_experiments]
    print("finish read", "time:", (time.time() - start_time) / 60, " m.")

    dist_info = []

    for i, experiment in enumerate(experiments):
        # Go through each step.
        for k in range(100, len(experiment), k_step):
            graph = experiment[k]
            # cur_dist_info = find_true_evolution_dist_fixed_parameters(
            #     graph, p_aa, p_bb, alpha, use_b_d_in_log=True
            # )
            cur_dist_info = find_true_evolution_dist_and_find_parameters1(
                graph, to_represent, use_b_d_in_log=True
            )

            cur_dist_info["experiment"] = i
            cur_dist_info["real_alpha"] = alpha
            cur_dist_info["real_p_aa"] = p_aa
            cur_dist_info["real_p_bb"] = p_bb
            cur_dist_info["empirical_true_dist"] = k
            dist_info.append(cur_dist_info)

        if i % 10 == 0:
            print("experiments:", i, "time:", (time.time() - start_time) / 60, " m.")

    log_dictionaries(
        dist_info,
        "3d_fragile_breakage_model/logs/" + method + file + ".csv",
    )

    print("time:", (time.time() - start_time) / 60, " m.")


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
            "plot": empirical_min_distances,
            "label": "Empirical min distance",
            "color": "black",
            "linestyle": "dashed",
        },
        {
            "plot": estimated_true_distances,
            "label": "Estimated true distance",
            "color": "blue",
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

    relative_error = []
    for i, true_distance in enumerate(empirical_true_distances):
        relative_error.append(
            100 * abs(true_distance - estimated_true_distances[i]) / true_distance
        )

    draw(
        xs,
        relative_error,
        "Real number of rearrangements",
        "% relative error",
        "Relative percentage error of estimated true evolutionary distance\n"
        + parameters_for_plot_title,
        save_as + "_relative_percentage_error",
    )


def draw_best_x(dist_info, parameters_for_plot_title, save_as):
    xs = list(
        map(
            lambda info: int(info["empirical_true_dist"])
            / parameters.NUMBER_OF_FRAGILE_EDGES,
            dist_info,
        )
    )
    best_xs = list(map(lambda info: float(info["best_x"]), dist_info))
    plots = [
        {"plot": xs, "label": "Empirical x", "color": "red", "linestyle": "dashed"},
        {
            "plot": best_xs,
            "label": "Estimated x",
            "color": "blue",
        },
    ]

    draw_plots(
        xs=xs,
        plots=plots,
        x_label="Real x",
        y_label="Estimated x",
        title="Estimated x\n" + parameters_for_plot_title,
        save_as=save_as,
    )


def draw_true_dist_for_parameters(folder_name, parameter_index, draw_parameters):
    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    file, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    p_ab = 1 - p_aa - p_bb
    plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
    print(plot_title)

    dist_info = read_logs(
        "3d_fragile_breakage_model/logs/" + folder_name + file + ".csv"
    )

    draw_dists(
        dist_info,
        None,
        None,
        plot_title,
        "3d_fragile_breakage_model/plots/" + folder_name + file,
    )

    draw_best_x(
        dist_info,
        plot_title,
        "3d_fragile_breakage_model/plots/" + folder_name + "best_x/" + file,
    )

    if draw_parameters:
        xs = list(
            map(
                lambda info: int(info["empirical_true_dist"])
                / parameters.NUMBER_OF_FRAGILE_EDGES,
                dist_info,
            )
        )
        plots = [
            {
                "plot": list(map(lambda info: float(info["real_alpha"]), dist_info)),
                "label": "Empirical alpha",
                "color": "red",
                "linestyle": "dashed",
            },
            {
                "plot": list(map(lambda info: float(info["best_alpha"]), dist_info)),
                "label": "Estimated alpha",
                "color": "blue",
            },
            {
                "plot": list(map(lambda info: float(info["real_p_aa"]), dist_info)),
                "label": "Empirical p_aa",
                "color": "purple",
                "linestyle": "dashed",
            },
            {
                "plot": list(map(lambda info: float(info["best_p_aa"]), dist_info)),
                "label": "Estimated p_aa",
                "color": "lightpink",
            },
            {
                "plot": [p_ab] * len(xs),
                "label": "Empirical p_ab",
                "color": "darkorange",
                "linestyle": "dashed",
            },
            {
                "plot": list(map(lambda info: float(info["best_p_ab"]), dist_info)),
                "label": "Estimated p_ab",
                "color": "dimgray",
            },
        ]

        draw_plots(
            xs=xs,
            plots=plots,
            x_label="Real x",
            y_label="Estimated parameters",
            title="Estimated parameters\n" + plot_title,
            save_as="3d_fragile_breakage_model/plots/"
            + folder_name
            + "best_parameters/"
            + file,
        )


def draw_true_dist_with_additional_plot(
    index_of_parameters, additional_folder, additional_label, folder_to_save
):
    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[index_of_parameters]
    parameters_str, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
    print(plot_title)

    all_dist_info_fixed = read_logs(
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
    need_ks = set(map(lambda info: info["empirical_true_dist"], dist_info_found))
    dist_info_fixed = []
    for k in range(len(all_dist_info_fixed)):
        if all_dist_info_fixed[k]["empirical_true_dist"] in need_ks:
            dist_info_fixed.append(all_dist_info_fixed[k])

    draw_dists(
        dist_info_fixed,
        dist_info_found,
        additional_label,
        plot_title,
        "3d_fragile_breakage_model/plots/" + folder_to_save + parameters_str,
    )


def get_median_errors_by_min_dist(parameter_index, number_of_experiments, k_step=100):
    max_cycle_len_with_types = 0
    start_time = time.time()
    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    file = parameter["parameters_str"]
    if parameter["experiments_in_one_bunch"] != 1:
        print("experiments_in_one_bunch should be 1")
        return

    print(file)

    experiments = read_experiments_cycles_info(
        get_experiments_dir(parameter["experiments_in_one_bunch"]) + file + ".csv",
        max_cycle_len_with_types,
        max_cycle_len_with_types,
        is_int=False,
        is_cycles_ordered=True,
    )[:number_of_experiments]
    print("finish read", "time:", (time.time() - start_time) / 60, " m.")

    median_errors_by_min_dist_k = []
    steps = len(experiments[0])

    for k in range(100, steps, k_step):
        sample = list(map(lambda experiment: (experiment[k].d - k) / k, experiments))
        median_errors_by_min_dist_k.append(
            {"k": k, "median_errors_by_min_dist": median(sample)}
        )

    log_dictionaries(
        median_errors_by_min_dist_k,
        "3d_fragile_breakage_model/logs/median_errors_by_min_dist/" + file + ".csv",
    )

    print(
        "finish get_median_errors_by_min_dist, time:",
        (time.time() - start_time) / 60,
        " m.",
    )


def draw_box_plot(folder_name, number_of_experiments, parameter_index, from_percentile, to_percentile):
    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    file, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
    print(plot_title)

    dist_info = read_logs(
        "3d_fragile_breakage_model/logs/" + folder_name + file + ".csv"
    )
    median_errors_by_min_dist_log = read_logs(
        "3d_fragile_breakage_model/logs/median_errors_by_min_dist/" + file + ".csv"
    )
    median_errors_by_min_dist_dict = {}
    for median_error in median_errors_by_min_dist_log:
        median_errors_by_min_dist_dict[median_error["k"]] = float(
            median_error["median_errors_by_min_dist"]
        )
    experiments = []
    cur_experiment_index = 0
    cur_experiment = []
    for experiment_x in dist_info:
        if int(experiment_x["experiment"]) != cur_experiment_index:
            cur_experiment_index = int(experiment_x["experiment"])
            experiments.append(cur_experiment)
            cur_experiment = [experiment_x]
        else:
            cur_experiment.append(experiment_x)
    if len(cur_experiment) > 0:
        experiments.append(cur_experiment)

    samples_k = []
    steps = len(experiments[0])
    ks = []
    median_errors_by_min_dist = []
    for k in range(steps):
        median_errors_by_min_dist.append(
            median_errors_by_min_dist_dict[experiments[0][k]["empirical_true_dist"]]
        )

        ks.append(
            float(experiments[0][k]["empirical_true_dist"])
            / parameters.NUMBER_OF_FRAGILE_EDGES
        )
        sample = list(
            map(
                lambda experiment: (
                    float(experiment[k]["estimated_true_dist"])
                    - float(experiment[k]["empirical_true_dist"])
                )
                / float(experiment[k]["empirical_true_dist"]),
                experiments,
            )
        )

        samples_k.append(sample)

    fig, ax1 = plt.subplots()
    ax1.boxplot(samples_k, sym="", whis=(from_percentile, to_percentile))
    ax1.set(
        axisbelow=True,
        title="Relative error in estimating evol. dist., whiskers: ["
        + str(from_percentile)
        + "%; "
        + str(to_percentile)
        + "%], boxes: 50%,\n"
        + plot_title
        + ", experiments = "
        + str(number_of_experiments),
        xlabel="x",
        ylabel="Relative error",
    )
    plt.grid()

    plt.scatter(
        range(1, steps + 1),
        median_errors_by_min_dist,
        label="Median relative error by min dist",
    )
    plt.legend(loc="lower left")

    ax1.set_xlim(0.5, len(samples_k) + 0.5)
    ax1.set_xticklabels(ks)

    save_as = "3d_fragile_breakage_model/plots/" + folder_name + file
    plt.savefig(save_as)
    plt.close()


def main():
    # method_to_run = "true_evolution_distance_found_parameters/1/log_"
    method_to_run = "true_evolution_distance_found_parameters/1/box_plot/log_"
    # method_to_run = "true_evolution_distance_fixed/box_plot/log_"
    index_of_parameters = 14
    number_of_experiments = 100
    k_step = 200
    # k_step = 100
    # compute_true_evolutionary_distance_for_box_plot(
    #     parameter_index=index_of_parameters,
    #     number_of_experiments=number_of_experiments,
    #     k_step=k_step,
    #     method=method_to_run,
    # )
    # get_median_errors_by_min_dist(index_of_parameters, number_of_experiments)
    from_percentile = 5
    to_percentile = 95
    draw_box_plot(
        folder_name=method_to_run,
        number_of_experiments=number_of_experiments,
        parameter_index=index_of_parameters,
        from_percentile=from_percentile,
        to_percentile=to_percentile
    )
    # index_of_parameters = 4
    # compute_true_evolutionary_distance(index_of_parameters, method_to_run)
    # draw_true_dist_for_parameters(
    #     method_to_run,
    #     parameter_index=index_of_parameters,
    #     draw_parameters=True,
    # )
    # draw_true_dist_with_additional_plot(
    #     index_of_parameters,
    #     "true_evolution_distance_found_parameters/1/log_",
    #     "Estimated true distance with found parameters",
    #     "true_evolution_distance_fixed_and_found_parameters/",
    # )


if __name__ == "__main__":
    main()
