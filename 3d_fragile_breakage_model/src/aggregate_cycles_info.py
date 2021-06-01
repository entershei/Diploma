import parameters
from compute_statistics import compute_d_by_cycles, compute_b_by_cycles
from utils import read_experiments_cycles_info, define_cycles_representative
from generate_directories_names import (
    create_new_directory_in_cycles_info,
    get_experiments_dir,
)
from utils import CyclesInfo, log_experiments


def sum_cycles_info(experiments, max_possible_cycles_len, cycles_representatives):
    # print("start sum")
    num_steps_of_markov_process = len(experiments[0])

    summed_cycles_info = []
    for i in range(num_steps_of_markov_process):
        sum_cnt_cycles_m = {}
        sum_cycles_with_edges_order = {}
        sum_a_in_non_trivial_cycles = 0
        sum_b_in_non_trivial_cycles = 0
        sum_a_in_non_trivial_cycles_part = 0
        sum_b_in_non_trivial_cycles_part = 0
        sum_d = 0
        sum_b = 0

        for cycle_type in cycles_representatives:
            sum_cycles_with_edges_order[cycle_type] = 0

        for c_len in range(1, max_possible_cycles_len):
            sum_cnt_cycles_m[str(c_len)] = 0

        for experiment in experiments:
            for cycle_type in cycles_representatives:
                if cycle_type in experiment[i].cycles_with_edges_order:
                    sum_cycles_with_edges_order[cycle_type] += experiment[
                        i
                    ].cycles_with_edges_order[cycle_type]
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
            sum_d += experiment[i].d
            sum_b += experiment[i].b

        summed_cycles_info.append(
            CyclesInfo(
                sum_cnt_cycles_m,
                sum_cycles_with_edges_order,
                sum_a_in_non_trivial_cycles,
                sum_b_in_non_trivial_cycles,
                sum_a_in_non_trivial_cycles_part,
                sum_b_in_non_trivial_cycles_part,
                sum_d,
                sum_b,
            )
        )

    # print("finish sum")
    return summed_cycles_info


# If a_in_non_trivial_cycles_part and b_in_non_trivial_cycles_part were not calculated for experiments, there will be
# negative values.
def aggregate_cycles_info(
    f_in,
    f_out,
    experiments_in_one_bunch,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    print("start read experiments")
    experiments = read_experiments_cycles_info(
        f_in,
        max_cycle_len_with_types,
        max_interesting_cycles_len,
        is_int=False,
        is_cycles_ordered=True,
    )
    print("finish read")

    _, cycles_representatives = define_cycles_representative(max_cycle_len_with_types)
    num_experiments = len(experiments) * experiments_in_one_bunch
    summed_cycles_info = sum_cycles_info(
        experiments, max_interesting_cycles_len, cycles_representatives
    )
    aggregated_cycles_info = []

    for info in summed_cycles_info:
        cycles_with_edges_order = {}
        for cycle_type in info.cycles_with_edges_order:
            cycles_with_edges_order[cycle_type] = (
                info.cycles_with_edges_order[cycle_type] / num_experiments
            )
        cycles_m = {}
        for cycle_len in info.cycles_m:
            cycles_m[cycle_len] = info.cycles_m[cycle_len] / num_experiments

        aggregated_cycles_info.append(
            CyclesInfo(
                cycles_m,
                cycles_with_edges_order,
                info.a_in_non_trivial_cycles / num_experiments,
                info.b_in_non_trivial_cycles / num_experiments,
                info.a_in_non_trivial_cycles_part / num_experiments,
                info.b_in_non_trivial_cycles_part / num_experiments,
                compute_d_by_cycles(cycles_m),
                compute_b_by_cycles(cycles_m),
            )
        )

    log_experiments(
        aggregated_cycles_info,
        f_out + ".csv",
        "w",
        max_interesting_cycles_len=max_interesting_cycles_len,
        cycles_representatives=cycles_representatives,
    )


def main():
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = parameters.MAX_POSSIBLE_CYCLES_LEN

    parameter_index = 6
    parameter = parameters.PROBABILITIES_WITH_ALPHA[parameter_index]

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
