import csv
import math

import matplotlib.pyplot as plt

import parameters
from utils import (
    generate_cycle_type,
    generate_cycle_types,
    create_new_directories_for_result_comparison,
    get_cycles_info_dir,
    get_relative_error_dir,
    get_analytical_cycles_dir,
)


def read_logs(f, cycles_types):
    num_c_n = []
    with open(f, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        i = 0
        for row in reader:
            n = int(row["n"])
            k = int(row["k"])
            assert i == k

            c_n = {}
            for cycle_type in cycles_types:
                c_n[cycle_type] = float(row[cycle_type])
            num_c_n.append(c_n)
            i += 1

    return num_c_n, n


# Возвращает число циклов длины m, посчитанных через аналитическую формулу.
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
            * alpha ** l
            * beta ** r
            * l ** (r - 1)
            / math.factorial(l)
            * r ** (l - 1)
            / math.factorial(r)
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
        / alpha ** (m - 2)
        / math.factorial(m)
        * math.exp(-x * m * (2 * p_aa + p_ab) / alpha)
    )

    all_b = (
        (2 * x * p_bb) ** (m - 1)
        * m ** (m - 2)
        / beta ** (m - 2)
        / math.factorial(m)
        * math.exp(-x * m * (2 * p_bb + p_ab) / beta)
    )

    cycles = {generate_cycle_type(m, m): all_a, generate_cycle_type(m, 0): all_b}

    sum_all = all_a + all_b
    for l in range(1, m):
        cur_cycles = cycles_depends_on_cnt_a(l)
        cycles[generate_cycle_type(m, l)] = cur_cycles
        sum_all += cur_cycles

    cycles["all"] = sum_all

    return cycles


# returns 100 * |analytical_c_n - real_c_n| / real_c_n
# alpha = t / n, t is number of A-type edges
def compute_relative_errors(
    cycle_len, x, real_c_ns, cycles_types, n, p_aa, p_bb, alpha
):
    k = int(x * n)

    relative_errors = {}
    real_c_n = {}
    all_real_c_n = 0
    for cycle_type in cycles_types:
        real_c_n[cycle_type] = real_c_ns[k][cycle_type] / n
        all_real_c_n += real_c_n[cycle_type]

    analytical_c_n = compute_analytical_cycles_m(cycle_len, x, p_aa, p_bb, alpha)

    for cycle_type in cycles_types:
        if real_c_n[cycle_type] != 0:
            relative_errors[cycle_type] = (
                100
                * abs(analytical_c_n[cycle_type] - real_c_n[cycle_type])
                / real_c_n[cycle_type]
            )
        else:
            relative_errors[cycle_type] = 0
    if all_real_c_n != 0:
        relative_errors["all"] = (
            100 * abs(analytical_c_n["all"] - all_real_c_n) / all_real_c_n
        )
    else:
        relative_errors["all"] = 0

    return relative_errors


def write_relative_error(error_depends_on_x, f, field_names):
    with open(f, "w", newline="") as log:
        fieldnames = ["x"] + field_names
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for x in error_depends_on_x.keys():
            error_depends_on_x[x]["x"] = x
            log_lens.writerow(error_depends_on_x[x])


def compare_m_cycles(cycle_len, f_in, f_out, p_aa, p_bb, alpha, cycles_types):
    real_c_n_s, n = read_logs(f_in, cycles_types)

    # x = k / n
    error_depends_on_x = {}
    x = 0.0
    step = 1 / n
    while x < 1 + step:
        error_depends_on_x[x] = compute_relative_errors(
            cycle_len,
            x,
            real_c_n_s,
            cycles_types,
            n,
            p_aa,
            p_bb,
            alpha,
        )
        x += step

    write_relative_error(error_depends_on_x, f_out, cycles_types + ["all"])


def draw_average(real_e, ess_e, xs):
    plt.plot(xs, real_e)
    plt.plot(xs, ess_e)
    plt.xlabel("Number of swaps")
    plt.ylabel("E")
    plt.grid()
    plt.show()


def write_analytical_cycles(analytical_cycles_depends_on_x, f, field_names):
    with open(f, "w", newline="") as log:
        fieldnames = field_names + ["x"]
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for x in analytical_cycles_depends_on_x.keys():
            analytical_cycles_depends_on_x[x]["x"] = x
            log_lens.writerow(analytical_cycles_depends_on_x[x])


def compute_analytical_cycles(n, cycle_len, p_aa, p_bb, alpha, f_out, field_names):
    # x = k / n
    analytical_cycles_depends_on_x = {}
    x = 0.0
    step = 1 / n
    while x < 1 + step:
        analytical_cycles_depends_on_x[x] = compute_analytical_cycles_m(
            cycle_len, x, p_aa, p_bb, alpha
        )
        x += step

    write_analytical_cycles(analytical_cycles_depends_on_x, f_out, field_names)


def result_comparison(cycle_len, file, p_aa, p_bb, alpha):
    file_end = file + parameters.EXPERIMENTS + ".csv"

    cycle_types = generate_cycle_types(cycle_len, cycle_len)

    compare_m_cycles(
        cycle_len=cycle_len,
        f_in=get_cycles_info_dir() + file_end,
        f_out=get_relative_error_dir(cycle_len) + "depends_on_x_" + file_end,
        p_aa=p_aa,
        p_bb=p_bb,
        alpha=alpha,
        cycles_types=cycle_types,
    )

    compute_analytical_cycles(
        n=parameters.NUMBER_OF_FRAGILE_EDGES,
        cycle_len=cycle_len,
        p_aa=p_aa,
        p_bb=p_bb,
        alpha=alpha,
        f_out=get_analytical_cycles_dir(cycle_len) + file + ".csv",
        field_names=cycle_types + ["all"],
    )


def main():
    max_interesting_cycles_len = 5
    create_new_directories_for_result_comparison(max_interesting_cycles_len)

    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        file, p_aa, p_bb, alpha = cur_parameters

        for i in range(max_interesting_cycles_len):
            result_comparison(
                i + 1,
                file,
                p_aa,
                p_bb,
                alpha,
            )


if __name__ == "__main__":
    main()
