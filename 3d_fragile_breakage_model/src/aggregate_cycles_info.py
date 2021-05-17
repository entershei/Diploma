import parameters
from utils import generate_cycle_types, read_experiments_cycles_info
from generate_directories_names import (
    create_new_directory_in_cycles_info,
    get_experiments_dir,
)
from utils import CyclesInfo, log_experiments


def sum_cycles_info(experiments, max_cycle_len_with_types, max_possible_cycles_len):
    # print("start sum")
    num_steps_of_markov_process = len(experiments[0])

    cycle_types = generate_cycle_types(1, max_cycle_len_with_types)

    summed_cycles_info = []
    for i in range(num_steps_of_markov_process):
        sum_cnt_cycle_types = {}
        sum_cnt_cycles_m = {}
        sum_a_in_non_trivial_cycles = 0
        sum_b_in_non_trivial_cycles = 0
        sum_a_in_non_trivial_cycles_part = 0
        sum_b_in_non_trivial_cycles_part = 0

        for cycle_type in cycle_types:
            sum_cnt_cycle_types[cycle_type] = 0

        for c_len in range(1, max_possible_cycles_len):
            sum_cnt_cycles_m[str(c_len)] = 0

        for experiment in experiments:
            for cycle_type in cycle_types:
                if cycle_type in experiment[i].cycle_types:
                    sum_cnt_cycle_types[cycle_type] += experiment[i].cycle_types[
                        cycle_type
                    ]
            for c_len in range(1, max_possible_cycles_len):
                if str(c_len) in experiment[i].cycles_m:
                    sum_cnt_cycles_m[str(c_len)] += experiment[i].cycles_m[str(c_len)]
            sum_a_in_non_trivial_cycles += experiment[i].a_in_non_trivial_cycles
            sum_b_in_non_trivial_cycles += experiment[i].b_in_non_trivial_cycles
            sum_a_in_non_trivial_cycles_part += experiment[
                i
            ].a_in_non_trivial_cycles_part
            sum_b_in_non_trivial_cycles_part += experiment[
                i
            ].b_in_non_trivial_cycles_part

        summed_cycles_info.append(
            CyclesInfo(
                sum_cnt_cycle_types,
                sum_cnt_cycles_m,
                sum_a_in_non_trivial_cycles,
                sum_b_in_non_trivial_cycles,
                sum_a_in_non_trivial_cycles_part,
                sum_b_in_non_trivial_cycles_part,
            )
        )

    # print("finish sum")
    return summed_cycles_info


# If a_in_non_trivial_cycles_part and b_in_non_trivial_cycles_part were not calculated for experiments, there will be
# negative values.
def aggregate_cycles_info(
    f_in, f_out, experiments_in_one_bunch, max_cycle_len_with_types, max_interesting_cycles_len
):
    print("start read experiments")
    experiments = read_experiments_cycles_info(
        f_in, max_cycle_len_with_types, max_interesting_cycles_len, is_int=True
    )
    print("finish read")

    num_experiments = len(experiments) * experiments_in_one_bunch
    summed_cycles_info = sum_cycles_info(
        experiments, max_cycle_len_with_types, max_interesting_cycles_len
    )
    aggregated_cycles_info = []

    for info in summed_cycles_info:
        cycle_types = {}
        for cycle_type in info.cycle_types:
            cycle_types[cycle_type] = info.cycle_types[cycle_type] / num_experiments
        cycles_m = {}
        for cycle_len in info.cycles_m:
            cycles_m[cycle_len] = info.cycles_m[cycle_len] / num_experiments

        aggregated_cycles_info.append(
            CyclesInfo(
                cycle_types,
                cycles_m,
                info.a_in_non_trivial_cycles / num_experiments,
                info.b_in_non_trivial_cycles / num_experiments,
                info.a_in_non_trivial_cycles_part / num_experiments,
                info.b_in_non_trivial_cycles_part / num_experiments,
            )
        )

    log_experiments(
        aggregated_cycles_info,
        f_out + ".csv",
        "w",
        max_cycle_len_with_types=max_cycle_len_with_types,
        max_interesting_cycles_len=max_interesting_cycles_len,
    )


def main():
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = parameters.MAX_POSSIBLE_CYCLES_LEN

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[-1:]:
        experiments_in_one_bunch = parameter["experiments_in_one_bunch"]
        cycles_info_log_path = create_new_directory_in_cycles_info(
            parameter["number_of_experiments"], experiments_in_one_bunch
        )

        string_parameters = parameter["parameters_str"]
        file = string_parameters + ".csv"

        print(string_parameters)

        aggregate_cycles_info(
            f_in=get_experiments_dir(experiments_in_one_bunch) + file,
            f_out=cycles_info_log_path + string_parameters,
            experiments_in_one_bunch=experiments_in_one_bunch,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )


if __name__ == "__main__":
    main()
