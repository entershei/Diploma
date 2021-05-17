from pathlib import Path
import parameters


def number_of_fragile_edges():
    return "n" + str(parameters.NUMBER_OF_FRAGILE_EDGES) + "/"


def get_parameters_as_string(number_of_experiments):
    return number_of_fragile_edges() + str(number_of_experiments) + "_experiments/"


def get_cycles_info_dir(number_of_experiments, experiments_in_one_bunch):
    return (
        "3d_fragile_breakage_model/logs/cycles_info/"
        + get_parameters_as_string(number_of_experiments)
        + "from_experiments_"
        + str(experiments_in_one_bunch)
    )


def create_new_directory_in_cycles_info(
    number_of_experiments, experiments_in_one_bunch
):
    path = get_cycles_info_dir(number_of_experiments, experiments_in_one_bunch)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def get_experiments_dir(experiments_in_one_bunch):
    return (
        "3d_fragile_breakage_model/logs/experiments_"
        + str(experiments_in_one_bunch)
        + "/"
        + number_of_fragile_edges()
    )


def create_new_directory_for_logging_experiments(experiments_in_one_bunch):
    path = get_experiments_dir(experiments_in_one_bunch)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path


def get_relative_error_dir(number_of_experiments, add_depends_on_k=True):
    path = "3d_fragile_breakage_model/logs/relative_error/" + get_parameters_as_string(
        number_of_experiments
    )
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


def get_plots_relative_error(number_of_experiments, folder_name):
    return (
        "3d_fragile_breakage_model/plots/relative_error/"
        + get_parameters_as_string(number_of_experiments)
        + folder_name
        + "/"
    )


def get_plots_compare_cycles(number_of_experiments, folder_name):
    return (
        "3d_fragile_breakage_model/plots/to_compare_number_of_cycles/"
        + get_parameters_as_string(number_of_experiments)
        + folder_name
        + "/"
    )


def create_new_directories_in_plots():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        folder = parameter["parameters_str"]
        number_of_experiments = parameter["number_of_experiments"]
        Path(get_plots_relative_error(number_of_experiments, folder)).mkdir(
            parents=True, exist_ok=True
        )
        Path(get_plots_compare_cycles(number_of_experiments, folder)).mkdir(
            parents=True, exist_ok=True
        )
