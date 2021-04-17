import csv

import matplotlib.pyplot as plt
import numpy as np

import parameters
from utils import (
    generate_cycle_types,
    create_new_directories_in_plots,
    get_plots_aggregated_cycles_dir,
    get_cycles_info_dir,
    get_plots_relative_error,
    get_plots_compare_cycles,
    get_relative_error_dir,
    get_analytical_cycles_dir,
)


def read_logs(f, cycle_types, to_sum, name_for_sum, to_rename):
    xs = []
    errors = []
    with open(f, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        for row in reader:
            xs.append(float(row["x"]))
            error = {}
            for cycle_type in cycle_types:
                error[cycle_type] = float(row[cycle_type])
            for r in to_rename:
                error[r[1]] = error[r[0]]

            summ = 0
            for cycle_type in to_sum:
                summ += error[cycle_type]
            error[name_for_sum] = summ

            errors.append(error)

    return xs, errors


def draw_error(xs, errors, title, save_as):
    plt.plot(xs, errors)
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("Relative error")
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def draw_relative_errors(
    folder_name,
    parameters_for_plot_name,
    number_of_experiments,
    max_cycle_len,
):
    for c_len in range(1, max_cycle_len + 1):
        cycle_len = str(c_len)
        cycle_names = generate_cycle_types(c_len, c_len)
        save_path = get_plots_relative_error(folder_name)
        xs, errors = read_logs(
            get_relative_error_dir(cycle_len)
            + "depends_on_x_"
            + folder_name
            + number_of_experiments
            + ".csv",
            cycle_names + ["all"],
            [],
            "",
            [["all", cycle_len]],
        )

        for cycle_type in cycle_names + [cycle_len]:
            draw_error(
                xs,
                list(map(lambda error: error[cycle_type], errors)),
                "Relative error of number of "
                + cycle_type
                + " depends on x,\n"
                + parameters_for_plot_name,
                save_path + cycle_type + number_of_experiments + ".png",
            )


def draw_number_of_cycles(cycles, title, save_as):
    k = len(cycles)
    x = list(range(k))
    plt.plot(x, cycles)
    plt.title(title)
    plt.xlabel("Number of swaps")
    plt.ylabel("Cycles")
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def draw_two_number_of_cycles(xs, real_cycles, analytical_cycles, title, save_as):
    plt.plot(xs, real_cycles, color="green", label="real")
    plt.plot(xs, analytical_cycles, color="blue", label="analytical")
    plt.legend(["Real", "Analytical"])
    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("Cycles")
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def number_of_cycles(cycle_len, arr):
    cycle_types = generate_cycle_types(cycle_len, cycle_len)
    res = np.add(arr[cycle_types[0]], arr[cycle_types[1]])
    for i in range(2, len(cycle_types)):
        res = np.add(res, arr[cycle_types[i]])
    return res


def read_cycles_info_logs(f, max_cycle_len):
    steps = 0
    num_all_cycles = []
    cycle_types = generate_cycle_types(1, 5)
    different_cycles = {}
    for cycle_type in cycle_types:
        different_cycles[cycle_type] = []

    with open(f, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        for row in reader:
            num_all_cycles.append(float(row["all"]))

            for cycle_type in cycle_types:
                different_cycles[cycle_type].append(float(row[cycle_type]))

            steps += 1

    cycles_info = {
        "num_all_cycles": num_all_cycles,
        "different_cycles": different_cycles,
        "steps": steps,
    }

    for i in range(1, max_cycle_len + 1):
        cycles_info[str(i)] = number_of_cycles(i, different_cycles)

    return cycles_info


def draw_average_cycles(folder_name, number_of_experiments, max_cycle_len):
    f = get_cycles_info_dir() + folder_name + number_of_experiments + ".csv"
    save_path = get_plots_aggregated_cycles_dir(folder_name)

    cycles_info = read_cycles_info_logs(f, max_cycle_len)

    draw_number_of_cycles(
        cycles_info["num_all_cycles"],
        "Average number of cycles depends of number of swaps",
        save_path + "all_cycles" + number_of_experiments + ".png",
    )

    for cycle_len in range(1, max_cycle_len + 1):
        cycle = str(cycle_len)
        draw_number_of_cycles(
            cycles_info[cycle],
            "Average number of " + cycle + "-cycles depends of number of swaps",
            save_path + cycle + "_cycles" + number_of_experiments + ".png",
        )

    cycle_types = generate_cycle_types(1, max_cycle_len)
    for cycle_type in cycle_types:
        draw_number_of_cycles(
            cycles_info["different_cycles"][cycle_type],
            "Average number of " + cycle_type + "-cycles depends of number of swaps",
            save_path
            + str(len(cycle_type))
            + "_"
            + cycle_type
            + "_cycles"
            + number_of_experiments
            + ".png",
        )


def interesting_cycles_info(n, xs, cycles_info):
    interesting_info = []
    for x in xs:
        interesting_info.append(cycles_info[int(x * n)] / n)
    return interesting_info


def draw_average_with_analytical_cycles(
    folder_name, number_of_experiments, max_cycle_len
):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    save_path = get_plots_compare_cycles(folder_name)

    f_real = get_cycles_info_dir() + folder_name + number_of_experiments + ".csv"

    cycles_info = read_cycles_info_logs(f_real, max_cycle_len)

    for c_len in range(1, max_cycle_len + 1):
        cycle_len = str(c_len)
        cycle_types = generate_cycle_types(c_len, c_len)
        xs, analytical_cycles = read_logs(
            get_analytical_cycles_dir(cycle_len) + folder_name + ".csv",
            cycle_types,
            cycle_types,
            cycle_len,
            [],
        )

        max_x = cycles_info["steps"] / n
        for i in range(len(xs)):
            if xs[i] > max_x:
                xs = xs[:i]
                analytical_cycles = analytical_cycles[:i]
                break

        draw_two_number_of_cycles(
            xs,
            interesting_cycles_info(n, xs, cycles_info[cycle_len]),
            list(map(lambda cycle: cycle[cycle_len], analytical_cycles)),
            "Normalized number of " + cycle_len + "-cycles depends of x",
            save_path + cycle_len + "_cycles" + number_of_experiments + ".png",
        )

        for cycle_type in cycle_types:
            draw_two_number_of_cycles(
                xs,
                interesting_cycles_info(
                    n, xs, cycles_info["different_cycles"][cycle_type]
                ),
                list(map(lambda cycle: cycle[cycle_type], analytical_cycles)),
                "Normalized number of " + cycle_type + "-cycles depends of x",
                save_path
                + cycle_len
                + "_"
                + cycle_type
                + "_cycles"
                + number_of_experiments
                + ".png",
            )


def main():
    max_interesting_cycles_len = 5
    create_new_directories_in_plots()

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        folder_name, p_aa, p_bb, alpha = cur_parameters
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

        draw_relative_errors(
            folder_name,
            parameters_for_plot_name,
            parameters.EXPERIMENTS,
            max_cycle_len=max_interesting_cycles_len,
        )
        draw_average_cycles(
            folder_name,
            parameters.EXPERIMENTS,
            max_cycle_len=max_interesting_cycles_len,
        )
        draw_average_with_analytical_cycles(
            folder_name,
            parameters.EXPERIMENTS,
            max_cycle_len=max_interesting_cycles_len,
        )


if __name__ == "__main__":
    main()
