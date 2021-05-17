import csv

import parameters
from aggregate_cycles_info import read_experiments_cycles_info
from compute_statistics import compute_analytical_cycles_m
from utils import (
    CyclesInfo,
    generate_cycle_types_for_len,
    log_dictionaries,
)
from generate_directories_names import (
    create_new_directories_for_result_comparison,
    get_cycles_info_dir,
    get_relative_error_dir,
    get_analytical_cycles_dir,
)


# Returns a_in_non_trivial_cycles, b_in_non_trivial_cycles only for for cycle_len = [2; max_interesting_cycles_len).
def compute_analytical_cycles(
    n,
    p_aa,
    p_bb,
    alpha,
    steps,
    f_out,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
):
    def write_analytical_cycles():
        with open(f_out, "w", newline="") as log:
            log_cycles = csv.DictWriter(log, fieldnames=field_names)
            log_cycles.writeheader()
            for cur_k in analytical_cycles_depends_on_k.keys():
                row = {
                    "k": cur_k,
                    "a_in_non_trivial_cycles": analytical_cycles_depends_on_k[
                        cur_k
                    ].a_in_non_trivial_cycles,
                    "b_in_non_trivial_cycles": analytical_cycles_depends_on_k[
                        cur_k
                    ].b_in_non_trivial_cycles,
                    "a_in_non_trivial_cycles_part": analytical_cycles_depends_on_k[
                        cur_k
                    ].a_in_non_trivial_cycles_part,
                    "b_in_non_trivial_cycles_part": analytical_cycles_depends_on_k[
                        cur_k
                    ].b_in_non_trivial_cycles_part,
                }
                for c_type in analytical_cycles_depends_on_k[cur_k].cycle_types.keys():
                    row[c_type] = analytical_cycles_depends_on_k[cur_k].cycle_types[
                        c_type
                    ]
                for c_len in analytical_cycles_depends_on_k[cur_k].cycles_m.keys():
                    row[c_len] = analytical_cycles_depends_on_k[cur_k].cycles_m[c_len]

                log_cycles.writerow(row)

    analytical_cycles_depends_on_k = {}
    field_names = [
        "k",
        "a_in_non_trivial_cycles",
        "b_in_non_trivial_cycles",
        "a_in_non_trivial_cycles_part",
        "b_in_non_trivial_cycles_part",
    ]
    k = 0
    while k < steps:
        x = k / n
        cycle_types = {}
        cycles_m = {}
        a_in_non_trivial_cycles = 0
        b_in_non_trivial_cycles = 0
        a_in_non_trivial_cycles_part = 0
        b_in_non_trivial_cycles_part = 0
        for cycle_len in range(1, max_interesting_cycles_len):
            cur_analytical_cycles = compute_analytical_cycles_m(
                cycle_len, x, p_aa, p_bb, alpha
            )
            if cycle_len < max_cycle_len_with_types:
                for cycle_type in cur_analytical_cycles["types"].keys():
                    cycle_types[cycle_type] = cur_analytical_cycles["types"][cycle_type]
                    if cycle_type not in field_names:
                        field_names.append(cycle_type)

            cycles_m[str(cycle_len)] = cur_analytical_cycles["all"]
            if str(cycle_len) not in field_names:
                field_names.append(str(cycle_len))

            if cycle_len > 1:
                for cycle_type in cur_analytical_cycles["types"].keys():
                    a_in_non_trivial_cycles = cur_analytical_cycles["types"][
                        cycle_type
                    ] * cycle_type.count("A")
                    b_in_non_trivial_cycles = cur_analytical_cycles["types"][
                        cycle_type
                    ] * cycle_type.count("B")

                    if cycle_len < parameters.PART:
                        a_in_non_trivial_cycles_part = cur_analytical_cycles["types"][
                            cycle_type
                        ] * cycle_type.count("A")
                        b_in_non_trivial_cycles_part = cur_analytical_cycles["types"][
                            cycle_type
                        ] * cycle_type.count("B")

        analytical_cycles_depends_on_k[k] = CyclesInfo(
            cycle_types,
            cycles_m,
            a_in_non_trivial_cycles,
            b_in_non_trivial_cycles,
            a_in_non_trivial_cycles_part,
            b_in_non_trivial_cycles_part,
        )
        k += 1

    write_analytical_cycles()
    return analytical_cycles_depends_on_k


# returns 100 * |analytical_cycles_info - empirical_cycles_info| / empirical_cycles_info
def compute_relative_error_between_two_results(
    empirical_cycles_info, analytical_cycles_info
):
    return (
        100
        * abs(empirical_cycles_info - analytical_cycles_info)
        / empirical_cycles_info
        if empirical_cycles_info != 0
        else 0
    )


def relative_error(
    empirical_cycles_info,
    analytical_cycles_info,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
    f_out,
):
    n = parameters.NUMBER_OF_FRAGILE_EDGES
    cycle_types = [[]]
    for i in range(1, max_cycle_len_with_types):
        cycle_types.append(generate_cycle_types_for_len(i))

    error_depends_on_k = []
    k = 0
    steps = len(empirical_cycles_info)
    while k < steps:
        errors = {"k": k}
        for cycle_len in range(1, max_interesting_cycles_len):
            if cycle_len < max_cycle_len_with_types:
                for cycle_type in cycle_types[cycle_len]:
                    errors[cycle_type] = compute_relative_error_between_two_results(
                        empirical_cycles_info[k].cycle_types[cycle_type] / n,
                        analytical_cycles_info[k].cycle_types[cycle_type],
                    )

            errors[str(cycle_len)] = compute_relative_error_between_two_results(
                empirical_cycles_info[k].cycles_m[str(cycle_len)] / n,
                analytical_cycles_info[k].cycles_m[str(cycle_len)],
            )
            errors[
                "a_in_non_trivial_cycles"
            ] = compute_relative_error_between_two_results(
                empirical_cycles_info[k].a_in_non_trivial_cycles / n,
                analytical_cycles_info[k].a_in_non_trivial_cycles,
            )
            errors[
                "b_in_non_trivial_cycles"
            ] = compute_relative_error_between_two_results(
                empirical_cycles_info[k].b_in_non_trivial_cycles / n,
                analytical_cycles_info[k].b_in_non_trivial_cycles,
            )
            if empirical_cycles_info[k].a_in_non_trivial_cycles_part >= 0:
                errors[
                    "a_in_non_trivial_cycles_part"
                ] = compute_relative_error_between_two_results(
                    empirical_cycles_info[k].a_in_non_trivial_cycles_part / n,
                    analytical_cycles_info[k].a_in_non_trivial_cycles_part,
                )
                errors[
                    "b_in_non_trivial_cycles_part"
                ] = compute_relative_error_between_two_results(
                    empirical_cycles_info[k].b_in_non_trivial_cycles_part / n,
                    analytical_cycles_info[k].b_in_non_trivial_cycles_part,
                )

        error_depends_on_k.append(errors)
        k += 1

    log_dictionaries(error_depends_on_k, f_out)


def main():
    max_interesting_cycles_len = 11
    max_cycle_len_with_types = 6

    create_new_directories_for_result_comparison()

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA[-1:]:
        file = cur_parameters["parameters_str"]
        number_of_experiments = cur_parameters["number_of_experiments"]

        empirical_cycles_info = read_experiments_cycles_info(
            get_cycles_info_dir(number_of_experiments, cur_parameters["experiments_in_one_bunch"]) + file + ".csv",
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            False,
        )[0]
        analytical_cycles_info = compute_analytical_cycles(
            n=parameters.NUMBER_OF_FRAGILE_EDGES,
            p_aa=cur_parameters["p_aa"],
            p_bb=cur_parameters["p_bb"],
            alpha=cur_parameters["alpha"],
            steps=len(empirical_cycles_info),
            f_out=get_analytical_cycles_dir() + file + ".csv",
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
        )

        relative_error(
            empirical_cycles_info=empirical_cycles_info,
            analytical_cycles_info=analytical_cycles_info,
            max_cycle_len_with_types=max_cycle_len_with_types,
            max_interesting_cycles_len=max_interesting_cycles_len,
            f_out=get_relative_error_dir(number_of_experiments) + file + ".csv",
        )


if __name__ == "__main__":
    main()
