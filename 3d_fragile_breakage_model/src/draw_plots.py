import matplotlib.pyplot as plt
import numpy as np
import parameters
from utils import (
    read_experiments_cycles_info,
    generate_cycle_types_for_len,
    define_cycles_representative,
)
from generate_directories_names import (
    create_new_directories_in_plots,
    get_cycles_info_dir,
    get_plots_relative_error,
    get_plots_compare_cycles,
    get_relative_error_dir,
    get_analytical_cycles_dir,
)


def draw(xs, ys, x_label, y_label, title, save_as):
    plt.plot(xs, ys)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid()
    plt.savefig(save_as)
    plt.close()


def draw_number_of_cycles(cycles, title, save_as):
    xs = list(range(len(cycles)))
    draw(xs, cycles, "Number of swaps", "Cycles", title, save_as)


def draw_plots(xs, plots, x_label, y_label, title, save_as):
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

    plt.legend()
    plt.title(title)
    plt.xlabel(x_label)
    if xs[-1] == 1.5:
        xs_positions = list(np.arange(0, max(xs) + 0.1, 0.2)) + [1.5]
        plt.xticks(xs_positions)
    plt.ylabel(y_label)
    plt.grid()
    plt.savefig(save_as)
    plt.close()


def draw_relative_errors(
    folder_name,
    parameters_for_plot_title,
    number_of_experiments,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
    to_represent,
):
    f_in_relative_error = (
        get_relative_error_dir(number_of_experiments) + folder_name + ".csv"
    )
    relative_errors = read_experiments_cycles_info(
        f_in_relative_error,
        max_cycle_len_with_types,
        max_interesting_cycles_len,
        is_int=False,
    )[0]

    n = parameters.NUMBER_OF_FRAGILE_EDGES
    ks = range(len(relative_errors))
    xs = list(map(lambda k: k / n, ks))

    save_path = get_plots_relative_error(number_of_experiments, folder_name)

    for c_len in range(1, max_interesting_cycles_len):
        cycle_len = str(c_len)

        if c_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(c_len):
                draw(
                    xs,
                    list(
                        map(
                            lambda k: relative_errors[k].cycles_with_edges_order[
                                to_represent[cycle_type]
                            ],
                            ks,
                        )
                    ),
                    "x",
                    "Relative error",
                    "Relative error of number of "
                    + to_represent[cycle_type]
                    + " depends on x,\n"
                    + parameters_for_plot_title,
                    save_path + cycle_type + ".png",
                )
        draw(
            xs,
            list(
                map(
                    lambda k: relative_errors[k].cycles_m[cycle_len],
                    ks,
                )
            ),
            "x",
            "Relative error",
            "Relative error of number of "
            + cycle_len
            + " depends on x,\n"
            + parameters_for_plot_title,
            save_path + cycle_len + ".png",
        )


def draw_cycles(
    get_empirical_cycles_lambda,
    get_analytical_cycles_lambda,
    xs,
    empirical_cycles_info,
    analytical_cycles_info,
    cycle_type,
    parameters_for_plot_title,
    save_path,
):
    draw_plots(
        xs,
        [
            {
                "plot": list(
                    map(
                        get_empirical_cycles_lambda,
                        empirical_cycles_info,
                    )
                ),
                "label": "Empirical cycles",
                "color": "red",
            },
            {
                "plot": list(
                    map(
                        get_analytical_cycles_lambda,
                        analytical_cycles_info,
                    )
                ),
                "label": "Analytical cycles",
                "color": "blue",
            },
        ],
        "x",
        "Cycles/n",
        title="Normalized number of "
        + cycle_type
        + "-cycles depends of x\n"
        + parameters_for_plot_title,
        save_as=save_path + cycle_type + "_cycles" + ".png",
    )


def draw_empirical_with_analytical_cycles(
    folder_name,
    parameters_for_plot_title,
    empirical_cycles_info,
    analytical_cycles_info,
    number_of_experiments,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
    to_represent,
):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    save_path = get_plots_compare_cycles(number_of_experiments, folder_name)
    xs = list(map(lambda k: k / n, range(len(empirical_cycles_info))))

    for c_len in range(1, max_interesting_cycles_len):
        if c_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(c_len):
                draw_cycles(
                    lambda cycle_info: cycle_info.cycles_with_edges_order[
                        to_represent[cycle_type]
                    ]
                    / n,
                    lambda cycle_info: cycle_info.cycles_with_edges_order[
                        to_represent[cycle_type]
                    ],
                    xs,
                    empirical_cycles_info,
                    analytical_cycles_info,
                    cycle_type,
                    parameters_for_plot_title,
                    save_path,
                )

        cycle_len = str(c_len)
        draw_cycles(
            lambda cycle_info: cycle_info.cycles_m[cycle_len] / n,
            lambda cycle_info: cycle_info.cycles_m[cycle_len],
            xs,
            empirical_cycles_info,
            analytical_cycles_info,
            cycle_len,
            parameters_for_plot_title,
            save_path,
        )


def build_parameters_for_plot_title(p_aa, p_bb, alpha):
    return (
        "n = "
        + str(parameters.NUMBER_OF_FRAGILE_EDGES)
        + ", p_aa = "
        + str(p_aa)
        + ", p_bb = "
        + str(p_bb)
        + ", α = "
        + str(alpha)
    )


def main():
    max_interesting_cycles_len = 11
    max_cycle_len_with_types = 6
    to_represent, _ = define_cycles_representative(max_cycle_len_with_types)
    create_new_directories_in_plots()

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA[-1:]:
        folder_name, p_aa, p_bb, alpha = (
            cur_parameters["parameters_str"],
            cur_parameters["p_aa"],
            cur_parameters["p_bb"],
            cur_parameters["alpha"],
        )
        print(folder_name)

        parameters_for_plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)

        f_in_empirical = (
            get_cycles_info_dir(
                cur_parameters["number_of_experiments"],
                cur_parameters["experiments_in_one_bunch"],
            )
            + folder_name
            + ".csv"
        )
        empirical_cycles_info = read_experiments_cycles_info(
            f_in_empirical,
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            is_int=False,
        )[0]

        f_in_analytically = get_analytical_cycles_dir() + folder_name + ".csv"
        analytical_cycles_info = read_experiments_cycles_info(
            f_in_analytically,
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            is_int=False,
        )[0]

        draw_relative_errors(
            folder_name,
            parameters_for_plot_title,
            cur_parameters["number_of_experiments"],
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
            to_represent=to_represent,
        )
        draw_empirical_with_analytical_cycles(
            folder_name,
            parameters_for_plot_title,
            empirical_cycles_info,
            analytical_cycles_info,
            cur_parameters["number_of_experiments"],
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
            to_represent=to_represent,
        )


def cycles_together(parameters_index):
    max_cycle_len_with_types = 1
    max_interesting_cycles_len = 7

    cur_parameters = parameters.PROBABILITIES_WITH_ALPHA[parameters_index]
    folder_name, p_aa, p_bb, alpha = (
        cur_parameters["parameters_str"],
        cur_parameters["p_aa"],
        cur_parameters["p_bb"],
        cur_parameters["alpha"],
    )
    number_of_experiments = cur_parameters["number_of_experiments"]
    print(folder_name)

    parameters_for_plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)

    f_in_empirical = (
        get_cycles_info_dir(
            number_of_experiments, cur_parameters["experiments_in_one_bunch"]
        )
        + folder_name
        + ".csv"
    )
    empirical_cycles_info = read_experiments_cycles_info(
        f_in_empirical,
        max_cycle_len_with_types,
        max_interesting_cycles_len,
        is_int=False,
    )[0]

    f_in_analytically = get_analytical_cycles_dir() + folder_name + ".csv"
    analytical_cycles_info = read_experiments_cycles_info(
        f_in_analytically,
        max_cycle_len_with_types,
        max_interesting_cycles_len,
        is_int=False,
    )[0]

    n = parameters.NUMBER_OF_FRAGILE_EDGES
    save_path = get_plots_compare_cycles(number_of_experiments, folder_name)
    max_steps = 1501
    xs = list(map(lambda k: k / n, range(len(empirical_cycles_info))))[:max_steps]
    empirical_cycles_info = empirical_cycles_info[:max_steps]
    analytical_cycles_info = analytical_cycles_info[:max_steps]

    draw_plots(
        xs,
        [
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["3"] / n,
                        empirical_cycles_info,
                    )
                ),
                "label": "Empirical c3/n",
                "color": "red",
            },
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["3"],
                        analytical_cycles_info,
                    )
                ),
                "label": "Analytical c3/n",
                "color": "blue",
                "linestyle": "dashed",
            },
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["4"] / n,
                        empirical_cycles_info,
                    )
                ),
                "label": "Empirical c4/n",
                "color": "dimgray",
            },
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["4"],
                        analytical_cycles_info,
                    )
                ),
                "label": "Analytical c4/n",
                "color": "lightpink",
                "linestyle": "dashed",
            },
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["5"] / n,
                        empirical_cycles_info,
                    )
                ),
                "label": "Empirical c5/n",
                "color": "darkorange",
            },
            {
                "plot": list(
                    map(
                        lambda cycle_info: cycle_info.cycles_m["5"],
                        analytical_cycles_info,
                    )
                ),
                "label": "Analytical c5/n",
                "color": "darkgreen",
                "linestyle": "dashed",
            },
        ],
        "x = steps/n",
        "Cycles/n",
        title="Normalized number of "
        + "cycles depends of x\n"
        + parameters_for_plot_title,
        save_as=save_path + "cycles_with_different_lens" + ".png",
    )


if __name__ == "__main__":
    # main()
    cycles_together(4)
