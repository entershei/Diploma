from pathlib import Path

import parameters


def generate_cycle_types(min_len, max_len):
    cycles = []
    for cur_len in range(min_len, max_len + 1):
        for cnt_a in range(cur_len, -1, -1):
            cycles.append("A" * cnt_a + "B" * (cur_len - cnt_a))
    return cycles


def number_of_fragile_edges():
    return "n" + str(parameters.NUMBER_OF_FRAGILE_EDGES) + "/"


def get_parameters_as_string():
    return number_of_fragile_edges() + "different_number_of_experiments/"


def get_cycles_info_dir():
    return "logs/cycles_info/" + get_parameters_as_string()


def create_new_directory_in_cycles_info():
    path = get_cycles_info_dir()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def create_new_directory_for_logging_experiments():
    path = "logs/experiments/" + number_of_fragile_edges()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def create_new_directories_for_result_comparison():
    for i in range(1, 4):
        Path(
            "logs/relative_error/" + get_parameters_as_string() + str(i) + "cycles"
        ).mkdir(parents=True, exist_ok=True)
        Path(
            "logs/analytical_cycles/" + number_of_fragile_edges() + "c" + str(i)
        ).mkdir(parents=True, exist_ok=True)


def create_new_directories_in_plots():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        folder = list(parameter)[0]
        Path("plots/aggregated_cycles/" + get_parameters_as_string() + folder).mkdir(
            parents=True, exist_ok=True
        )
        Path("plots/relative_error/" + get_parameters_as_string() + folder).mkdir(
            parents=True, exist_ok=True
        )
        Path(
            "plots/to_compare_number_of_cycles/" + get_parameters_as_string() + folder
        ).mkdir(parents=True, exist_ok=True)
