import csv
import math

import parameters
from aggregate_cycles_info import read_experiments_cycles_info
from utils import (
    CyclesInfo,
    generate_cycle_type,
    generate_cycle_types_for_len,
    log_dictionaries,
)
from generate_directories_names import (
    create_new_directories_for_result_comparison,
    get_cycles_info_dir,
    get_relative_error_dir,
    get_analytical_cycles_dir,
)


# Возвращает нормированное число циклов длины m, посчитанных через аналитическую формулу, указывает число циклов по
# типам (например, AAAB-циклы) и их общее количество.
def compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha):
    p_ab = 1 - p_aa - p_bb
    beta = 1 - alpha

    def cycles_depends_on_cnt_a(l):
        eps_zero = 1e-9
        if p_ab < eps_zero or alpha < eps_zero or beta < eps_zero:
            return 0.0

        r = m - l

        return (
            (x * p_ab / (alpha * beta)) ** (m - 1)
            * (alpha * r) ** (l - 1)
            / math.factorial(r)
            * (beta * l) ** (r - 1)
            / math.factorial(l)
            * alpha
            * beta
            * (1 + 2 * beta * l * p_aa / (alpha * r * p_ab)) ** (l - 1)
            * (1 + 2 * alpha * r * p_bb / (beta * l * p_ab)) ** (r - 1)
            * math.exp(
                -x
                * (
                    2 * beta * l * p_aa
                    + alpha * r * p_ab
                    + beta * l * p_ab
                    + 2 * alpha * r * p_bb
                )
                / (alpha * beta)
            )
        )

    all_a = (
        (2 * x * p_aa) ** (m - 1)
        * m ** (m - 2)
        / math.factorial(m)
        / alpha ** (m - 2)
        * math.exp(-x * m * (2 * p_aa + p_ab) / alpha)
    )

    all_b = (
        (2 * x * p_bb) ** (m - 1)
        * m ** (m - 2)
        / math.factorial(m)
        / beta ** (m - 2)
        * math.exp(-x * m * (2 * p_bb + p_ab) / beta)
    )

    cycles = {generate_cycle_type(m, m): all_a, generate_cycle_type(m, 0): all_b}

    sum_all = all_a + all_b
    for l in range(1, m):
        cur_cycles = cycles_depends_on_cnt_a(l)
        cycles[generate_cycle_type(m, l)] = cur_cycles
        sum_all += cur_cycles

    return {"types": cycles, "all": sum_all}


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

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        file, p_aa, p_bb, alpha = cur_parameters

        empirical_cycles_info = read_experiments_cycles_info(
            get_cycles_info_dir() + file + ".csv",
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            False,
        )[0]
        analytical_cycles_info = compute_analytical_cycles(
            n=parameters.NUMBER_OF_FRAGILE_EDGES,
            p_aa=p_aa,
            p_bb=p_bb,
            alpha=alpha,
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
            f_out=get_relative_error_dir() + file + ".csv",
        )


if __name__ == "__main__":
    main()
