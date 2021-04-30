import csv

import matplotlib.pyplot as plt

import parameters
from utils import (
    generate_cycle_types,
    read_experiments_cycles_info,
    parse_logs_row,
    generate_cycle_types_for_len,
)
from generate_directories_names import (
    create_new_directories_in_plots,
    get_plots_aggregated_cycles_dir,
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
        plt.plot(xs, plot["plot"], label=plot["label"], color=plot["color"])
    plt.legend()
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid()
    plt.savefig(save_as)
    plt.close()


def read_logs_with_x(f, max_cycle_len_with_types, max_interesting_cycles_len):
    possible_cycle_types = generate_cycle_types(1, max_cycle_len_with_types)
    experiment = {}

    with open(f, "r", newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        for row in reader:
            x = float(row["x"])
            cycle_info = parse_logs_row(
                row, possible_cycle_types, max_interesting_cycles_len, False
            )
            experiment[x] = cycle_info

        csvfile.close()

    return experiment


def draw_relative_errors(
    folder_name,
    parameters_for_plot_title,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    f_in_relative_error = get_relative_error_dir() + folder_name + ".csv"
    relative_errors = read_logs_with_x(
        f_in_relative_error, max_cycle_len_with_types, max_interesting_cycles_len
    )

    save_path = get_plots_relative_error(folder_name)

    for c_len in range(1, max_interesting_cycles_len):
        cycle_len = str(c_len)

        if c_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(c_len):
                draw(
                    relative_errors.keys(),
                    list(
                        map(
                            lambda x: relative_errors[x].cycle_types[cycle_type],
                            relative_errors.keys(),
                        )
                    ),
                    "x",
                    "Relative error",
                    "Relative error of number of "
                    + cycle_type
                    + " depends on x,\n"
                    + parameters_for_plot_title,
                    save_path + cycle_type + ".png",
                )
        draw(
            relative_errors.keys(),
            list(
                map(
                    lambda x: relative_errors[x].cycles_m[cycle_len],
                    relative_errors.keys(),
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


def draw_average_cycles(
    folder_name,
    parameters_for_plot_title,
    cycles_info,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    save_path = get_plots_aggregated_cycles_dir(folder_name)

    for cycle_len in range(1, max_interesting_cycles_len):
        if cycle_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(cycle_len):
                draw_number_of_cycles(
                    list(
                        map(
                            lambda cycle_info: cycle_info.cycle_types[cycle_type],
                            cycles_info,
                        )
                    ),
                    "Average number of "
                    + cycle_type
                    + "-cycles depends of number of swaps\n"
                    + parameters_for_plot_title,
                    save_path + cycle_type + "_cycles" + ".png",
                )
        draw_number_of_cycles(
            list(
                map(
                    lambda cycle_info: cycle_info.cycles_m[str(cycle_len)],
                    cycles_info,
                )
            ),
            "Average number of "
            + str(cycle_len)
            + "-cycles depends of number of swaps\n"
            + parameters_for_plot_title,
            save_path + str(cycle_len) + "_cycles" + ".png",
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
                        analytical_cycles_info.keys(),
                    )
                ),
                "label": "Analytical cycles",
                "color": "blue",
            },
        ],
        "x",
        "Cycles",
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
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    save_path = get_plots_compare_cycles(folder_name)
    xs = analytical_cycles_info.keys()

    for c_len in range(1, max_interesting_cycles_len):
        if c_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(c_len):
                draw_cycles(
                    lambda cycle_info: cycle_info.cycle_types[cycle_type] / n,
                    lambda x: analytical_cycles_info[x].cycle_types[cycle_type],
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
            lambda x: analytical_cycles_info[x].cycles_m[cycle_len],
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
        + ", paa = "
        + str(p_aa)
        + ", p_bb = "
        + str(p_bb)
        + ", α = "
        + str(alpha)
    )


def main():
    max_interesting_cycles_len = 11
    max_cycle_len_with_types = 6

    create_new_directories_in_plots()

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        folder_name, p_aa, p_bb, alpha = cur_parameters
        print(folder_name)

        parameters_for_plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)

        f_in_empirical = get_cycles_info_dir() + folder_name + ".csv"
        empirical_cycles_info = read_experiments_cycles_info(
            f_in_empirical, max_cycle_len_with_types, max_interesting_cycles_len, False
        )[0]

        f_in_analytically = get_analytical_cycles_dir() + folder_name + ".csv"
        analytical_cycles_info = read_logs_with_x(
            f_in_analytically, max_cycle_len_with_types, max_interesting_cycles_len
        )

        draw_relative_errors(
            folder_name,
            parameters_for_plot_title,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )
        draw_average_cycles(
            folder_name,
            parameters_for_plot_title,
            empirical_cycles_info,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )
        draw_empirical_with_analytical_cycles(
            folder_name,
            parameters_for_plot_title,
            empirical_cycles_info,
            analytical_cycles_info,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )


if __name__ == "__main__":
    main()
