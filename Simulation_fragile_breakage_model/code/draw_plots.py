import csv

import matplotlib.pyplot as plt
import numpy as np

import parameters
from utils import generate_cycle_names, generate_cycle_types


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
    file_end, experiments, parameters_for_plot_name, max_cycle_len
):
    for c_len in range(1, max_cycle_len + 1):
        cycle_len = str(c_len)
        cycle_names = generate_cycle_names(c_len, c_len)
        save_path = "plots/relative_error/" + experiments + file_end + "/"
        xs, errors = read_logs(
            "logs/relative_error/"
            + experiments
            + cycle_len
            + "cycles/"
            + "depends_on_x_"
            + file_end
            + ".csv",
            cycle_names + ["all"],
            [],
            "",
            [["all", cycle_len + "-cycles"]],
        )

        for cycle_type in cycle_names + [cycle_len + "-cycles"]:
            draw_error(
                xs,
                list(map(lambda error: error[cycle_type], errors)),
                "Relative error of number of "
                + cycle_type
                + " depends on x,\n"
                + parameters_for_plot_name,
                save_path + cycle_type + ".png",
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


def read_cycles_info_logs(f):
    num_all_cycles = []
    cycle_types = generate_cycle_types(1, 5)
    different_cycles = {}
    for cycle_type in cycle_types:
        different_cycles[cycle_type] = []

    with open(f, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        for row in reader:
            num_all_cycles.append(float(row["all-cycles"]))

            for cycle_type in cycle_types:
                different_cycles[cycle_type].append(float(row[cycle_type + "-cycles"]))

    return {
        "num_all_cycles": num_all_cycles,
        "different_cycles": different_cycles,
        "1-cycles": number_of_cycles(1, different_cycles),
        "2-cycles": number_of_cycles(2, different_cycles),
        "3-cycles": number_of_cycles(3, different_cycles),
    }


def draw_average_cycles(file_end, experiments, max_cycle_len):
    parameters_string = experiments + file_end
    f = "logs/cycles_info/" + parameters_string + ".csv"
    save_path = "plots/aggregated_cycles/" + parameters_string + "/"

    cycles_info = read_cycles_info_logs(f)

    draw_number_of_cycles(
        cycles_info["num_all_cycles"],
        "Average number of cycles depends of number of swaps",
        save_path + "all_cycles.png",
    )

    for cycle_len in range(1, max_cycle_len + 1):
        cycle = str(cycle_len)
        draw_number_of_cycles(
            cycles_info[cycle + "-cycles"],
            "Average number of " + cycle + "-cycles depends of number of swaps",
            save_path + cycle + "_cycles.png",
        )

    cycle_types = generate_cycle_types(1, 3)
    for cycle_type in cycle_types:
        draw_number_of_cycles(
            cycles_info["different_cycles"][cycle_type],
            "Average number of " + cycle_type + "-cycles depends of number of swaps",
            save_path + str(len(cycle_type)) + "_" + cycle_type + "_cycles.png",
        )


def interesting_cycles_info(n, xs, cycles_info):
    interesting_info = []
    for x in xs:
        interesting_info.append(cycles_info[int(x * n)] / n)
    return interesting_info


def draw_average_with_analytical_cycles(file_end, experiments, max_cycle_len):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    real_parameters = experiments + file_end
    save_path = "plots/to_compare_number_of_cycles/" + real_parameters + "/"

    f_real = "logs/cycles_info/" + real_parameters + ".csv"

    cycles_info = read_cycles_info_logs(f_real)

    for c_len in range(1, max_cycle_len + 1):
        cycle_len = str(c_len)
        cycle_types = generate_cycle_names(cycle_len, cycle_len)
        xs, analytical_cycles = read_logs(
            "logs/analytical_cycles/n"
            + str(parameters.NUMBER_OF_FRAGILE_EDGES)
            + "/c"
            + cycle_len
            + "/"
            + file_end
            + ".csv",
            cycle_types,
            cycle_types,
            cycle_len + "-cycles",
            [],
        )
        draw_two_number_of_cycles(
            xs,
            interesting_cycles_info(n, xs, cycles_info[cycle_len + "-cycles"]),
            list(map(lambda cycle: cycle[cycle_len + "-cycles"], analytical_cycles)),
            "Normalized number of " + cycle_len + "-cycles depends of x",
            save_path + cycle_len + "_cycles.png",
        )

        for cycle_type in cycle_types:
            draw_two_number_of_cycles(
                xs,
                interesting_cycles_info(
                    n, xs, cycles_info["different_cycles"][cycle_type]
                ),
                list(
                    map(lambda cycle: cycle[cycle_type + "-cycles"], analytical_cycles)
                ),
                "Normalized number of " + cycle_type + "-cycles depends of x",
                save_path + cycle_len + "_" + cycle_type + "_cycles.png",
            )


def main():
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        file_end, parameters_for_plot_name = cur_parameters
        draw_relative_errors(
            file_end,
            parameters.EXPERIMENTS_FIELD_NAME,
            parameters_for_plot_name,
            max_cycle_len=3,
        )
        draw_average_cycles(
            file_end, parameters.EXPERIMENTS_FIELD_NAME, max_cycle_len=3
        )
        draw_average_with_analytical_cycles(
            file_end, parameters.EXPERIMENTS_FIELD_NAME, max_cycle_len=3
        )


if __name__ == "__main__":
    main()
