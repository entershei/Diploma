import time

import parameters
from generate_directories_names import get_cycles_info_dir
from compute_statistics import (
    compute_analytically_d_n,
    compute_analytically_b_n,
    estimate_alpha,
    estimate_p_ab,
    compute_empirical_d,
    compute_empirical_b,
)
from draw_plots import draw_plots, build_parameters_for_plot_title
from src.compute_statistics import compute_analytical_cycles_m
from utils import read_experiments_cycles_info, log_dictionaries, read_logs
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


def find_true_evolution_dist_fixed_parameters(graph, p_aa, p_bb, alpha):
    empirical_d = compute_empirical_d(graph)
    empirical_b = compute_empirical_b(graph)

    step = 5e-4
    l = step
    r = 1.5 + step

    min_error = 1e10
    best_x = 0.0

    for x in np.arange(l, r, step):
        error = compute_error(x, p_aa, p_bb, alpha, graph, empirical_b)
        if error < min_error:
            best_x = x
            min_error = error

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


def find_true_evolution_dist_with_parameters_using_skopt(graph, n_calls):
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
            graph.cycle_types["AA"],
            graph.cycle_types["AB"],
            graph.cycle_types["BB"],
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

    best_alphas = estimate_alpha(
        best_x,
        best_p_aa,
        best_p_bb,
        graph.cycle_types["AA"],
        graph.cycle_types["AB"],
        graph.cycle_types["BB"],
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


def find_true_evolution_dist_and_find_parameters0(graph):
    empirical_d = compute_empirical_d(graph)
    empirical_b = compute_empirical_b(graph)

    step_x = 1e-3
    l_x = step_x
    r_x = 1.5 + step_x

    min_error = 1e10

    best_analytical_x = 0.0
    best_p_aa = 0.0
    best_p_ab = 0.0
    best_p_bb = 0.0
    best_alpha = 0.0
    p_step = 0.015

    for x in np.arange(l_x, r_x, step_x):
        for p_aa in np.arange(p_step, 1, p_step):
            for alpha in np.arange(p_step, 1, p_step):
                p_ab = estimate_p_ab(
                    x, p_aa, alpha, graph.cycle_types["AA"], graph.cycle_types["AAA"]
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


def find_true_evolution_dist_and_find_parameters1(graph):
    empirical_d = compute_empirical_d(graph)
    empirical_b = compute_empirical_b(graph)

    step_x = 5e-4
    l_x = step_x
    r_x = 1.5 + step_x

    min_error = 1e10

    best_x = 0.0
    best_p_aa = 0.0
    best_p_ab = 0.0
    best_p_bb = 0.0
    best_alpha = 0.0
    p_step = 0.01

    for x in np.arange(l_x, r_x, step_x):
        for p_aa in np.arange(p_step, 1, p_step):
            for p_ab in np.arange(0, min(1 / 3, (1 - p_aa) / 2), p_step):
                p_bb = max(1 - p_aa - p_ab, 0.0)
                alphas = estimate_alpha(
                    x,
                    p_aa,
                    p_bb,
                    graph.cycle_types["AA"],
                    graph.cycle_types["AB"],
                    graph.cycle_types["BB"],
                )

                for alpha in alphas:
                    error = compute_error(x, p_aa, p_bb, alpha, graph, empirical_b)
                    if error < min_error:
                        best_x = x
                        best_p_aa = p_aa
                        best_p_ab = p_ab
                        best_p_bb = p_bb
                        best_alpha = alpha
                        min_error = error

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
        get_cycles_info_dir(parameter["number_of_experiments"]) + file + ".csv",
        max_cycle_len_with_types,
        parameters.MAX_POSSIBLE_CYCLES_LEN,
        False,
    )[0]

    dist_info = []

    # Go through each step.
    for k, graph in enumerate(graphs):
        if k < 1 or (k % 100) != 0:
            continue

        # cur_dist_info = find_true_evolution_dist_different_errors(graph, p_aa, p_bb, alpha)
        # cur_dist_info = find_true_evolution_dist_fixed_parameters(
        #     graph, p_aa, p_bb, alpha
        # )

        cur_dist_info = find_true_evolution_dist_and_find_parameters0(graph)
        # cur_dist_info = find_true_evolution_dist_and_find_parameters1(graph)
        # cur_dist_info = find_true_evolution_dist_with_parameters_using_skopt(
        #     graph, 2000
        # )
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
                (time.time() - start_time) / 60 / 60,
                " h.",
            )

    log_dictionaries(
        dist_info,
        "3d_fragile_breakage_model/logs/true_evolution_distance_found_parameters/"
        + method
        + "/"
        + file
        + ".csv",
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
            "plot": list(map(lambda v: float(v["estimated_true_dist"]), to_draw_main)),
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


def draw_true_dist_for_parameters(f_name, parameter_index, draw_parameters):
    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    folder_name, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    p_ab = 1 - p_aa - p_bb
    plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
    print(plot_title)

    dist_info = read_logs(
        "3d_fragile_breakage_model/logs/" + f_name + folder_name + ".csv"
    )

    draw_dists(
        dist_info,
        None,
        None,
        plot_title,
        "3d_fragile_breakage_model/plots/" + f_name + folder_name,
    )

    draw_best_x(
        dist_info,
        plot_title,
        "3d_fragile_breakage_model/plots/" + f_name + "best_x/" + folder_name,
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
            # {
            #     "plot": list(map(lambda info: float(info["real_p_bb"]), dist_info)),
            #     "label": "Empirical p_bb",
            #     "color": "darkorange",
            #     "linestyle": "dashed",
            # },
            # {
            #     "plot": list(map(lambda info: float(info["best_p_bb"]), dist_info)),
            #     "label": "Estimated p_bb",
            #     "color": "dimgray",
            # },
        ]

        draw_plots(
            xs=xs,
            plots=plots,
            x_label="Real x",
            y_label="Estimated parameters",
            title="Estimated parameters\n" + plot_title,
            save_as="3d_fragile_breakage_model/plots/"
            + f_name
            + "best_parameters/"
            + folder_name,
        )


def draw_true_dist_with_additional_plot(
    additional_folder, additional_label, folder_to_save
):
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        parameters_str, p_aa, p_bb, alpha = (
            cur_parameters["parameters_str"],
            cur_parameters["p_aa"],
            cur_parameters["p_bb"],
            cur_parameters["alpha"],
        )
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
    compute_true_evolutionary_distance(parameter_index=4, method="0")
    draw_true_dist_for_parameters(
        "true_evolution_distance_found_parameters/0/",
        parameter_index=4,
        draw_parameters=True,
    )
    # draw_true_dist_with_additional_plot(
    #     "true_evolution_distance_found_parameters/",
    #     "Estimated true distance with found parameters",
    #     "true_evolution_distance_fixed_and_found_parameters/",
    # )
