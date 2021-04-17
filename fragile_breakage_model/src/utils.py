from pathlib import Path

import parameters


def generate_cycle_type(cycle_len, cnt_a):
    return "A" * cnt_a + "B" * (cycle_len - cnt_a)


def generate_cycle_types(min_len, max_len):
    cycles = []
    for cur_len in range(min_len, max_len + 1):
        for cnt_a in range(cur_len, -1, -1):
            cycles.append(generate_cycle_type(cur_len, cnt_a))
    return cycles


def number_of_fragile_edges():
    return "n" + str(parameters.NUMBER_OF_FRAGILE_EDGES) + "/"


def get_parameters_as_string():
    return number_of_fragile_edges() + "different_number_of_experiments/"


def get_cycles_info_dir():
    return "fragile_breakage_model/logs/cycles_info/" + get_parameters_as_string()


def create_new_directory_in_cycles_info():
    path = get_cycles_info_dir()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def create_new_directory_for_logging_experiments():
    path = "fragile_breakage_model/logs/experiments/" + number_of_fragile_edges()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def get_relative_error_dir(cycle_len):
    return (
        "fragile_breakage_model/logs/relative_error/"
        + get_parameters_as_string()
        + str(cycle_len)
        + "cycles/"
    )


def get_analytical_cycles_dir(cycle_len):
    return (
        "fragile_breakage_model/logs/analytical_cycles/"
        + number_of_fragile_edges()
        + "c"
        + str(cycle_len)
        + "/"
    )


def create_new_directories_for_result_comparison():
    for i in range(1, 4):
        Path(get_relative_error_dir(i)).mkdir(parents=True, exist_ok=True)
        Path(get_analytical_cycles_dir(i)).mkdir(parents=True, exist_ok=True)


def create_new_directories_in_plots():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        folder = list(parameter)[0]
        Path(
            "fragile_breakage_model/plots/aggregated_cycles/"
            + get_parameters_as_string()
            + folder
        ).mkdir(parents=True, exist_ok=True)
        Path(
            "fragile_breakage_model/plots/relative_error/"
            + get_parameters_as_string()
            + folder
        ).mkdir(parents=True, exist_ok=True)
        Path(
            "fragile_breakage_model/plots/to_compare_number_of_cycles/"
            + get_parameters_as_string()
            + folder
        ).mkdir(parents=True, exist_ok=True)
