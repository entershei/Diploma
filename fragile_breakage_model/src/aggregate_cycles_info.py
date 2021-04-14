import parameters
import csv
from utils import generate_cycle_types, create_new_directory_in_cycles_info
from main import CyclesInfo


def read_experiments(f_in, max_cycle_len=5):
    cycle_types = generate_cycle_types(1, max_cycle_len)
    experiments = []

    with open(f_in, "r", newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        experiment = []
        for row in reader:
            k = int(row["k"])
            all_cycles = int(row["all"])
            max_cycle_len = int(row["max_cycle_len"])
            c_n = {}
            for cycle_type in cycle_types:
                c_n[cycle_type] = int(row[cycle_type])

            cycle_info = CyclesInfo(all_cycles, c_n, max_cycle_len)

            if k == 0 and len(experiment) > 0:
                if len(experiments) % 500 == 0:
                    print(k, len(experiments))
                experiments.append(experiment)
                experiment = [cycle_info]
            else:
                experiment.append(cycle_info)

        experiments.append(experiment)

        csvfile.close()

    return experiments


def aggregate_cycles_info(f_in, f_out, max_cycle_len=5):
    print("start read experiments")
    experiments = read_experiments(f_in, max_cycle_len)
    print("finish read")

    num_experiments = len(experiments)
    num_steps_of_markov_process = len(experiments[0])

    aggregated_cycles_info = []
    for i in range(num_steps_of_markov_process):
        average_all_cycles_num = 0
        average_cnt_cycles = {}
        for cycle_type in generate_cycle_types(1, max_cycle_len):
            average_cnt_cycles[cycle_type] = 0

        average_max_cycles_len = 0
        for experiment in experiments:
            average_all_cycles_num += experiment[i].num_all_cycles
            for cycles_types in experiment[i].num_n_cycles.keys():
                average_cnt_cycles[cycles_types] += experiment[i].num_n_cycles[
                    cycles_types
                ]

            average_max_cycles_len += experiment[i].max_len

        average_all_cycles_num /= num_experiments
        for cycles_types in average_cnt_cycles.keys():
            average_cnt_cycles[cycles_types] /= num_experiments
        average_max_cycles_len /= num_experiments

        aggregated_cycles_info.append(
            CyclesInfo(
                average_all_cycles_num, average_cnt_cycles, average_max_cycles_len
            )
        )

    print("start_log")

    log_aggregated_results(
        parameters.NUMBER_OF_FRAGILE_EDGES,
        aggregated_cycles_info,
        f_out + "_" + str(num_experiments) + "_experiments.csv",
    )
    print("finish log")


def log_aggregated_results(n, aggregated_cycles_info, file, max_cycle_len=5):
    with open(file, "w", newline="") as f_log_lens:
        cycle_types = generate_cycle_types(1, max_cycle_len)
        fieldnames = ["n", "k", "all"] + cycle_types + ["max_cycle_len"]

        log_cycles_info = csv.DictWriter(f_log_lens, fieldnames=fieldnames)
        log_cycles_info.writeheader()

        for step, info in enumerate(aggregated_cycles_info):
            cur_result = {
                "n": n,
                "k": step,
                "all": info.num_all_cycles,
                "max_cycle_len": info.max_len,
            }
            for cycle_type in cycle_types:
                cur_result[cycle_type] = info.num_n_cycles[cycle_type]
            log_cycles_info.writerow(cur_result)


def main():
    cycles_info_log_path = create_new_directory_in_cycles_info()

    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        string_parameters, p_aa, p_bb, a_type_edges_proportion = parameter
        file = string_parameters + ".csv"

        print(string_parameters)

        experiments_folder = (
            "logs/experiments/n" + str(parameters.NUMBER_OF_FRAGILE_EDGES) + "/"
        )

        aggregate_cycles_info(
            f_in=experiments_folder + file,
            f_out=cycles_info_log_path + string_parameters,
        )


if __name__ == "__main__":
    main()
