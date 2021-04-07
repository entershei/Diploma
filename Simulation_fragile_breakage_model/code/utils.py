from pathlib import Path

import parameters


def generate_cycle_types(min_len, max_len):
    cycles = []
    for cur_len in range(min_len, max_len + 1):
        for cnt_a in range(cur_len, -1, -1):
            cycles.append("A" * cnt_a + "B" * (cur_len - cnt_a))
    return cycles


def get_parameters_as_string():
    return (
        "n"
        + str(parameters.NUMBER_OF_FRAGILE_EDGES)
        + "/"
        + str(parameters.DIFFERENT_FRAGILE_EDGES_SPLITS)
        + "_"
        + str(parameters.MARKOV_PROCESS_EXPERIMENTS_ON_ONE_SPLIT)
        + "_experiments"
    )


def create_new_directory_in_cycles_info():
    Path("logs/cycles_info/" + get_parameters_as_string()).mkdir(
        parents=True, exist_ok=True
    )


def create_new_directories_in_relative_error_logs():
    for i in range(1, 4):
        Path(
            "logs/relative_error/"
            + get_parameters_as_string()
            + "/"
            + str(i)
            + "cycles"
        ).mkdir(parents=True, exist_ok=True)


def create_new_directories_in_plots():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        folder = list(parameter)[0][:-4]
        Path(
            "plots/aggregated_cycles/" + get_parameters_as_string() + "/" + folder
        ).mkdir(parents=True, exist_ok=True)
        Path("plots/relative_error/" + get_parameters_as_string() + "/" + folder).mkdir(
            parents=True, exist_ok=True
        )
        Path(
            "plots/to_compare_number_of_cycles/"
            + get_parameters_as_string()
            + "/"
            + folder
        ).mkdir(parents=True, exist_ok=True)
