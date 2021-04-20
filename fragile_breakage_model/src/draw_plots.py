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


def draw_two_plots(xs, ys1, ys2, x_label, y_label, title, legend, save_as):
    plt.plot(xs, ys1, color="green", label=legend[0])
    plt.plot(xs, ys2, color="blue", label=legend[1])
    plt.legend(legend)
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
    parameters_for_plot_name,
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
                    + parameters_for_plot_name,
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
            + parameters_for_plot_name,
            save_path + cycle_len + ".png",
        )


def draw_average_cycles(
    folder_name,
    parameters_for_plot_name,
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
                    + parameters_for_plot_name,
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
            + parameters_for_plot_name,
            save_path + str(cycle_len) + "_cycles" + ".png",
        )


def draw_empirical_with_analytical_cycles(
    folder_name,
    parameters_for_plot_name,
    empirical_cycles_info,
    analytical_cycles_info,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    save_path = get_plots_compare_cycles(folder_name)
    xs = list(
        map(
            lambda x: int(x * n),
            analytical_cycles_info.keys(),
        )
    )

    for c_len in range(1, max_interesting_cycles_len):
        if c_len < max_cycle_len_with_types:
            for cycle_type in generate_cycle_types_for_len(c_len):
                draw_two_plots(
                    xs,
                    list(
                        map(
                            lambda cycle_info: cycle_info.cycle_types[cycle_type] / n,
                            empirical_cycles_info,
                        )
                    ),
                    list(
                        map(
                            lambda x: analytical_cycles_info[x].cycle_types[cycle_type],
                            analytical_cycles_info.keys(),
                        )
                    ),
                    "x",
                    "Cycles",
                    "Normalized number of "
                    + cycle_type
                    + "-cycles depends of x\n"
                    + parameters_for_plot_name,
                    ["Real", "Analytical"],
                    save_path + cycle_type + "_cycles" + ".png",
                )

        cycle_len = str(c_len)
        draw_two_plots(
            xs,
            list(
                map(
                    lambda cycle_info: cycle_info.cycles_m[cycle_len] / n,
                    empirical_cycles_info,
                )
            ),
            list(
                map(
                    lambda x: analytical_cycles_info[x].cycles_m[cycle_len],
                    analytical_cycles_info.keys(),
                )
            ),
            "x",
            "Cycles",
            "Normalized number of "
            + cycle_len
            + "-cycles depends of x\n"
            + parameters_for_plot_name,
            ["Real", "Analytical"],
            save_path + cycle_len + "_cycles" + ".png",
        )


def main():
    max_interesting_cycles_len = 11
    max_cycle_len_with_types = 6

    create_new_directories_in_plots()

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA[:5]:
        folder_name, p_aa, p_bb, alpha = cur_parameters
        print(folder_name)

        parameters_for_plot_name = (
            "n = "
            + str(parameters.NUMBER_OF_FRAGILE_EDGES)
            + ", paa = "
            + str(p_aa)
            + ", p_bb = "
            + str(p_bb)
            + ", α = "
            + str(alpha)
        )

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
            parameters_for_plot_name,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )
        draw_average_cycles(
            folder_name,
            parameters_for_plot_name,
            empirical_cycles_info,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )
        draw_empirical_with_analytical_cycles(
            folder_name,
            parameters_for_plot_name,
            empirical_cycles_info,
            analytical_cycles_info,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )


if __name__ == "__main__":
    main()
