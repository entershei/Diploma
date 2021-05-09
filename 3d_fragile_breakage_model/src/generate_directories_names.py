from pathlib import Path
import parameters


def number_of_fragile_edges():
    return "n" + str(parameters.NUMBER_OF_FRAGILE_EDGES) + "/"


def get_parameters_as_string():
    return (
        number_of_fragile_edges()
        + str(parameters.NUMBER_OF_EXPERIMENTS)
        + "_experiments/"
    )


def get_cycles_info_dir():
    return "3d_fragile_breakage_model/logs/cycles_info/" + get_parameters_as_string()


def create_new_directory_in_cycles_info():
    path = get_cycles_info_dir()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def get_experiments_dir():
    return (
        "3d_fragile_breakage_model/logs/experiments_"
        + str(parameters.EXPERIMENTS_IN_ONE_BUNCH)
        + "/"
        + number_of_fragile_edges()
    )


def create_new_directory_for_logging_experiments():
    path = get_experiments_dir()
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def get_relative_error_dir(add_depends_on_k=True):
    path = "3d_fragile_breakage_model/logs/relative_error/" + get_parameters_as_string()
    if add_depends_on_k:
        return path + "depends_on_k_"
    else:
        return path


def get_analytical_cycles_dir():
    return (
        "3d_fragile_breakage_model/logs/analytical_cycles/" + number_of_fragile_edges()
    )


def create_new_directories_for_result_comparison():
    Path(get_relative_error_dir(False)).mkdir(parents=True, exist_ok=True)
    Path(get_analytical_cycles_dir()).mkdir(parents=True, exist_ok=True)


def get_plots_relative_error(folder_name):
    return (
        "3d_fragile_breakage_model/plots/relative_error/"
        + get_parameters_as_string()
        + folder_name
        + "/"
    )


def get_plots_compare_cycles(folder_name):
    return (
        "3d_fragile_breakage_model/plots/to_compare_number_of_cycles/"
        + get_parameters_as_string()
        + folder_name
        + "/"
    )


def create_new_directories_in_plots():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        folder = list(parameter)[0]
        Path(get_plots_relative_error(folder)).mkdir(parents=True, exist_ok=True)
        Path(get_plots_compare_cycles(folder)).mkdir(parents=True, exist_ok=True)
