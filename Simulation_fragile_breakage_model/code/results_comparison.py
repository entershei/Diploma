import csv
import math

import matplotlib.pyplot as plt

import parameters
from utils import generate_cycle_types


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


def compute_analytical_c1(x, p_aa, p_bb, alpha):
    beta = 1 - alpha
    a_cycles = alpha * math.exp(-1 * x * (1 + p_aa - p_bb) / alpha)
    b_cycles = beta * math.exp(-1 * x * (1 + p_bb - p_aa) / beta)
    return {"A": a_cycles, "B": b_cycles, "all": a_cycles + b_cycles}


def compute_analytical_c2(x, p_aa, p_bb, alpha):
    beta = 1 - alpha
    p_ab = 1 - p_aa - p_bb
    aa_cycles = p_aa * x * math.exp(-2 * x * (1 + p_aa - p_bb) / alpha)
    bb_cycles = p_bb * x * math.exp(-2 * x * (1 + p_bb - p_aa) / beta)
    ab_cycles = (
        p_ab
        * x
        * math.exp(
            -1
            * x
            * (2 * alpha * p_bb + alpha * p_ab + 2 * beta * p_aa + beta * p_ab)
            / (alpha * beta)
        )
    )

    return {
        "AA": aa_cycles,
        "AB": ab_cycles,
        "BB": bb_cycles,
        "all": aa_cycles + ab_cycles + bb_cycles,
    }


def compute_analytical_c3(x, p_aa, p_bb, alpha):
    beta = 1 - alpha
    p_ab = 1 - p_aa - p_bb

    aaa_cycles = (
        2 * (p_aa ** 2) * x ** 2 / alpha * math.exp(-3 * x * (2 * p_aa + p_ab) / alpha)
    )
    aab_cycles = (
        2 * x ** 2 * p_aa * p_ab / alpha + (p_ab * x) ** 2 / (2 * beta)
    ) * math.exp(
        -1
        * x
        * (2 * alpha * p_bb + alpha * p_ab + 4 * beta * p_aa + 2 * beta * p_ab)
        / (alpha * beta)
    )
    abb_cycles = (
        2 * x ** 2 * p_bb * p_ab / beta + (p_ab * x) ** 2 / (2 * alpha)
    ) * math.exp(
        -1
        * x
        * (2 * beta * p_aa + beta * p_ab + 4 * alpha * p_bb + 2 * alpha * p_ab)
        / (alpha * beta)
    )
    bbb_cycles = (
        2 * (p_bb ** 2) * x ** 2 / beta * math.exp(-3 * x * (2 * p_bb + p_ab) / beta)
    )

    return {
        "AAA": aaa_cycles,
        "AAB": aab_cycles,
        "ABB": abb_cycles,
        "BBB": bbb_cycles,
        "all": aaa_cycles + aab_cycles + abb_cycles + bbb_cycles,
    }


# returns 100 * |analytical_c_n - real_c_n| / real_c_n
# alpha = t / n, t is number of A-type edges
def compute_relative_errors(
    x, real_c_ns, cycles_types, n, p_aa, p_bb, alpha, estimation_func
):
    k = int(x * n)

    relative_errors = {}
    real_c_n = {}
    all_real_c_n = 0
    for cycle_type in cycles_types:
        real_c_n[cycle_type] = real_c_ns[k][cycle_type] / n
        all_real_c_n += real_c_n[cycle_type]

    analytical_c_n = estimation_func(x, p_aa, p_bb, alpha)

    for cycle_type in cycles_types:
        if real_c_n[cycle_type] != 0:
            relative_errors[cycle_type] = (
                100
                * abs(analytical_c_n[cycle_type] - real_c_n[cycle_type])
                / real_c_n[cycle_type]
            )
        else:
            relative_errors[cycle_type] = 0
    relative_errors["all"] = (
        100 * abs(analytical_c_n["all"] - all_real_c_n) / all_real_c_n
    )

    return relative_errors


def write_relative_error(error_depends_on_x, f, field_names):
    with open(f, "w", newline="") as log:
        fieldnames = ["x"] + field_names
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for x in error_depends_on_x.keys():
            error_depends_on_x[x]["x"] = x
            log_lens.writerow(error_depends_on_x[x])


def compare_n_cycles(
    f_in, f_out, p_aa, p_bb, alpha, cycles_types, compute_analytical_c_n
):
    real_c_n_s, n = read_logs(f_in, cycles_types)

    # x = k / n
    error_depends_on_x = {}
    x = 0.05
    step = 0.05
    while x < 2.05:
        error_depends_on_x[x] = compute_relative_errors(
            x, real_c_n_s, cycles_types, n, p_aa, p_bb, alpha, compute_analytical_c_n
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


def compute_analytical_cycles(
    p_aa, p_bb, alpha, compute_analytical_c_n, f_out, field_names
):
    # x = k / n
    analytical_cycles_depends_on_x = {}
    x = 0.05
    step = 0.05
    while x < 2.05:
        analytical_cycles_depends_on_x[x] = compute_analytical_c_n(x, p_aa, p_bb, alpha)

        x += step

    write_analytical_cycles(analytical_cycles_depends_on_x, f_out, field_names)


def result_comparison(
    cycles, file_end, p_aa, p_bb, alpha, experiments, compute_analytical_c_n
):
    cycle_types = generate_cycle_types(cycles, cycles)

    compare_n_cycles(
        f_in="logs/cycles_info/" + experiments + file_end,
        f_out="logs/relative_error/"
        + experiments
        + str(cycles)
        + "cycles/depends_on_x_"
        + file_end,
        p_aa=p_aa,
        p_bb=p_bb,
        alpha=alpha,
        cycles_types=cycle_types,
        compute_analytical_c_n=compute_analytical_c_n,
    )

    compute_analytical_cycles(
        p_aa=p_aa,
        p_bb=p_bb,
        alpha=alpha,
        compute_analytical_c_n=compute_analytical_c_n,
        f_out="logs/analytical_cycles/n"
        + str(parameters.NUMBER_OF_FRAGILE_EDGES)
        + "/c"
        + str(cycles)
        + "/"
        + file_end,
        field_names=cycle_types + ["all"],
    )


def main():
    for cur_parameters in parameters.PROBABILITIES_WITH_ALPHA:
        file_end, p_aa, p_bb, alpha = cur_parameters

        result_comparison(
            1,
            file_end,
            p_aa,
            p_bb,
            alpha,
            parameters.EXPERIMENTS_FIELD_NAME,
            compute_analytical_c1,
        )
        result_comparison(
            2,
            file_end,
            p_aa,
            p_bb,
            alpha,
            parameters.EXPERIMENTS_FIELD_NAME,
            compute_analytical_c2,
        )
        result_comparison(
            3,
            file_end,
            p_aa,
            p_bb,
            alpha,
            parameters.EXPERIMENTS_FIELD_NAME,
            compute_analytical_c3,
        )


if __name__ == "__main__":
    main()
