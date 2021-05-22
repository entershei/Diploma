import time

from compute_statistics import (
    compute_analytical_cycles_m,
    compute_empirical_b,
    compute_analytically_b_n,
)
from draw_plots import build_parameters_for_plot_title
from generate_directories_names import (
    get_experiments_dir,
    create_directory_for_confidence_intervals,
    get_confidence_intervals_dir,
    get_confidence_intervals_dir_plots,
)
import parameters
from utils import read_experiments_cycles_info, log_dictionaries, read_logs
from statistics import mean
import matplotlib.pyplot as plt
import numpy as np


def get_intervals(
    parameter_index,
    cycle_len,
    left_quantile,
    right_quantile,
    divide,
    number_of_experiments,
):
    start_time = time.time()

    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    experiments_in_one_bunch = parameter["experiments_in_one_bunch"]
    if experiments_in_one_bunch != 1:
        print("Wrong number of experiments in one bunch")
        exit(0)

    string_parameters = parameter["parameters_str"]

    print(string_parameters)

    f_in = get_experiments_dir(experiments_in_one_bunch) + string_parameters + ".csv"

    max_interesting_cycles_len = (
        parameters.MAX_POSSIBLE_CYCLES_LEN if divide == "b" else cycle_len + 1
    )
    experiments = read_experiments_cycles_info(
        f_in,
        0,
        max_interesting_cycles_len,
        is_int=True,
        number_of_experiments=number_of_experiments,
        is_cycles_ordered=True,
    )
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    steps = len(experiments[0])
    intervals_info = []

    for k in range(steps):
        if k == 0:
            continue
        sample_x = list(
            map(
                lambda experiment: experiment[k].cycles_m[str(cycle_len)]
                / (n if divide == "n" else compute_empirical_b(experiment[k])),
                experiments,
            )
        )

        sample_x.sort()

        interval_begins_index = int(len(sample_x) * left_quantile)
        intervals_ends_index = len(sample_x) - int(len(sample_x) * right_quantile) - 1
        x = k / n
        p_aa, p_bb, alpha = parameter["p_aa"], parameter["p_bb"], parameter["alpha"]
        intervals_info.append(
            {
                "x": x,
                "interval_begins": sample_x[interval_begins_index],
                "interval_ends": sample_x[intervals_ends_index],
                "mean": mean(sample_x),
                "analytical": compute_analytical_cycles_m(
                    cycle_len,
                    x,
                    p_aa,
                    p_bb,
                    alpha,
                )["all"]
                / (
                    1
                    if divide == "n"
                    else compute_analytically_b_n(x, p_aa, p_bb, alpha)
                ),
            }
        )

        if k % 100 == 0:
            print(k, (time.time() - start_time) / 60, " m.")

    directory = create_directory_for_confidence_intervals(
        len(experiments), cycle_len, divide
    )
    log_dictionaries(
        intervals_info,
        directory
        + build_file_name(string_parameters, left_quantile, right_quantile)
        + ".csv",
    )


def build_file_name(string_parameters, left_quantile, right_quantile):
    return (
        string_parameters
        + "_p:"
        + convert_quantile_to_str(left_quantile)
        + "_"
        + convert_quantile_to_str(right_quantile)
    )


def convert_quantile_to_str(quantile):
    return str(quantile * 100).replace(".", ",")


def draw_many_confidence_interval(
    parameter_index,
    cycles_len,
    left_quantile,
    right_quantile,
    colors,
    divide,
    number_of_experiments,
    between=False,
):
    def cycles_to_str():
        res = ""
        for j, c in enumerate(cycles_len):
            res += str(c)
            if j != len(cycles_len) - 1:
                res += "_"
        return res

    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]
    string_parameters = parameter["parameters_str"]
    plots = []
    xs = []
    intervals = []
    for i, cycle_len in enumerate(cycles_len):
        cur_intervals_info = read_logs(
            get_confidence_intervals_dir(number_of_experiments, cycle_len, divide)
            + build_file_name(string_parameters, left_quantile, right_quantile)
            + ".csv"
        )
        intervals.append(cur_intervals_info)
        intervals_begins = list(
            map(lambda info: float(info["interval_begins"]), cur_intervals_info)
        )
        intervals_ends = list(
            map(lambda info: float(info["interval_ends"]), cur_intervals_info)
        )
        plots.append(
            {
                "plot": list(map(lambda info: float(info["mean"]), cur_intervals_info)),
                "label": "Mean empirical с" + str(cycle_len) + "/" + divide,
                "color": colors[i * 3],
            }
        )
        plots.append(
            {
                "plot": list(
                    map(lambda info: float(info["analytical"]), cur_intervals_info)
                ),
                "label": "Analytical с" + str(cycle_len) + "/" + divide,
                "color": colors[i * 3 + 1],
                "linestyle": "dashed",
            },
        )

        xs = list(map(lambda info: float(info["x"]), cur_intervals_info))
        plt.fill_between(
            xs,
            intervals_begins,
            intervals_ends,
            color=colors[i * 3 + 2],
            label="c" + str(cycle_len) + "/" + divide + " confidence interval",
        )

    if between:
        between_xs = []
        cur_interval_begins = []
        cur_interval_ends = []
        first_begins = list(
            map(lambda info: float(info["interval_begins"]), intervals[0])
        )
        second_ends = list(map(lambda info: float(info["interval_ends"]), intervals[1]))
        for i in range(len(xs)):
            if first_begins[i] < second_ends[i]:
                between_xs.append(xs[i])
                cur_interval_begins.append(second_ends[i])
                cur_interval_ends.append(first_begins[i])

        plt.fill_between(
            between_xs,
            cur_interval_begins,
            cur_interval_ends,
            color=colors[-1],
        )

    for i, plot in enumerate(plots):
        if "linestyle" not in plot:
            plot["linestyle"] = "solid"
        plt.plot(
            xs,
            plot["plot"],
            label=plot["label"],
            color=plot["color"],
            linestyle=plot["linestyle"],
        )

    save_as = get_confidence_intervals_dir_plots(
        number_of_experiments, cycles_to_str(), divide
    ) + build_file_name(string_parameters, left_quantile, right_quantile)

    plt.legend()
    plt.title(
        "Normalized number of "
        + "m-cycles, "
        + "confidence interval from "
        + convert_quantile_to_str(left_quantile)
        + "% to "
        + convert_quantile_to_str(1 - right_quantile)
        + "%\n"
        + build_parameters_for_plot_title(
            parameter["p_aa"], parameter["p_bb"], parameter["alpha"]
        )
    )
    plt.xlabel("x")
    if xs[-1] == 1.5:
        xs_positions = list(np.arange(0, max(xs) + 0.1, 0.2)) + [1.5]
        plt.xticks(xs_positions)
    plt.ylabel("Cycles/" + divide)
    plt.grid()
    plt.savefig(save_as)
    plt.close()


def build_confidence_intervals(
    parameter_index,
    cycle_len,
    left_quantile,
    right_quantile,
    divide,
    number_of_experiments,
):
    get_intervals(
        parameter_index,
        cycle_len,
        left_quantile,
        right_quantile,
        divide,
        number_of_experiments,
    )
    draw_many_confidence_interval(
        parameter_index,
        [cycle_len],
        left_quantile,
        right_quantile,
        [
            "darkviolet",
            "blue",
            "lightblue",
        ],
        divide,
        number_of_experiments,
    )


def main():
    left_quantile = right_quantile = 0.05
    divide = "n"
    number_of_experiments = 1000
    index = 14

    build_confidence_intervals(
        index, 2, left_quantile, right_quantile, divide, number_of_experiments
    )
    build_confidence_intervals(
        index, 3, left_quantile, right_quantile, divide, number_of_experiments
    )
    build_confidence_intervals(
        index, 4, left_quantile, right_quantile, divide, number_of_experiments
    )
    build_confidence_intervals(
        index, 5, left_quantile, right_quantile, divide, number_of_experiments
    )
    build_confidence_intervals(
        index, 6, left_quantile, right_quantile, divide, number_of_experiments
    )

    colors = [
        "darkviolet",
        "blue",
        "lightblue",
        "orange",
        "darkgreen",
        "moccasin",
        "maroon",
        "lightslategray",
        "lightgray",
        "#b7e1a1",
        # "lightgray",
    ]
    draw_many_confidence_interval(
        index,
        [2, 4],
        left_quantile,
        right_quantile,
        colors,
        divide,
        number_of_experiments,
        between=True,
    )
    draw_many_confidence_interval(
        index,
        [3, 6],
        left_quantile,
        right_quantile,
        colors,
        divide,
        number_of_experiments,
        between=True,
    )


if __name__ == "__main__":
    main()
