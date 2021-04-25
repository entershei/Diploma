from results_comparison import compute_analytical_cycles_m
from generate_directories_names import get_cycles_info_dir
import matplotlib.pyplot as plt

import parameters
from aggregate_cycles_info import read_experiments_cycles_info
from draw_plots import draw
import numpy as np


# Normalized min distance
def compute_analytically_d_n(x, p_aa, p_bb, alpha, max_m):
    dist = 1.0

    for m in range(1, max_m):
        dist -= compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]

    return dist


# For non trivial cycles normalized sum of cycles lens
def compute_analytically_b_n(x, p_aa, p_bb, alpha):
    return 1 - compute_analytical_cycles_m(1, x, p_aa, p_bb, alpha)["all"]


#                 "real_true_dist": k,
#                 "analytical_true_dist": analytical_true_dist,
#                 "real_min_dist": empirical_min_d,
def draw_dists(to_draw, save_as):
    xs = list(map(lambda v: v["real_true_dist"], to_draw))
    empirical_true_distances = xs
    empirical_min_distances = list(map(lambda v: v["real_min_dist"], to_draw))
    analytical_true_distances = list(map(lambda v: v["analytical_true_dist"], to_draw))

    legend = [
        "Empirical true distance",
        "Estimated true distance",
        "Empirical min distance",
    ]

    plt.plot(xs, empirical_true_distances, color="blue", label=legend[0])
    plt.plot(xs, analytical_true_distances, color="red", label=legend[1])
    plt.plot(xs, empirical_min_distances, color="black", label=legend[2])

    plt.title("Evolutionary distance")
    plt.xlabel("Real number of rearrangements")
    plt.ylabel("Estimated number of rearrangements")
    plt.legend(legend)
    plt.grid()
    plt.savefig(save_as)
    plt.close()


def find_true_evolution_dist(parameters_str, p_aa, p_bb, alpha, max_m):
    def d_divide_by_b(possible_x):
        analytical_min_d = compute_analytically_d_n(
            possible_x,
            p_aa,
            p_bb,
            alpha,
            max_m,
        )
        analytical_sum_lens = compute_analytically_b_n(possible_x, p_aa, p_bb, alpha)
        return analytical_min_d / analytical_sum_lens

    aggregated_info = read_experiments_cycles_info(
        get_cycles_info_dir() + parameters_str + ".csv",
        6,
        parameters.MAX_POSSIBLE_CYCLES_LEN,
        False,
    )[0]

    print(parameters_str)

    step = 5e-4
    l = step
    r = 1 + step

    to_draw = []
    for k, info in enumerate(aggregated_info):
        if k == 0 or k % 5 != 0:
            continue
        print(k, "!")

        empirical_d = 0
        empirical_b = 0
        for cycle_len in info.cycles_m.keys():
            if int(cycle_len) > 1:
                empirical_d += info.cycles_m[cycle_len] * (int(cycle_len) - 1)
                empirical_b += info.cycles_m[cycle_len] * int(cycle_len)

        empirical_d_b = empirical_d / empirical_b

        best_analytical_x = 0
        min_error = 10000
        for x in np.arange(l, r, step):
            analytical_d_b = d_divide_by_b(x)
            cur_error = abs(analytical_d_b - empirical_d_b)
            if cur_error < min_error:
                best_analytical_x, min_error = x, cur_error
        analytical_b = compute_analytically_b_n(best_analytical_x, p_aa, p_bb, alpha)
        analytical_n = empirical_b / analytical_b
        analytical_true_dist = best_analytical_x * analytical_n

        to_draw.append(
            {
                "real_true_dist": k,
                "analytical_true_dist": analytical_true_dist,
                "real_min_dist": empirical_d,
            }
        )

    draw_dists(
        to_draw,
        "fragile_breakage_model/plots/true_evolution_distance/" + parameters_str,
    )


def compute_d_divide_b(x, p_aa, p_bb, alpha, max_m):
    analytical_min_d = compute_analytically_d_n(x, p_aa, p_bb, alpha, max_m)
    b = compute_analytically_b_n(x, p_aa, p_bb, alpha)
    if b > 0:
        return analytical_min_d / b
    return 0


def check_monotone_d_divide_b():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        s, p_aa, p_bb, a_type_edges_proportion = parameter
        print(s)

        max_m = 85
        d_b_s = []

        for k in range(1, parameters.NUMBER_OF_FRAGILE_EDGES):
            x = k / parameters.NUMBER_OF_FRAGILE_EDGES
            d_b = compute_d_divide_b(x, p_aa, p_bb, a_type_edges_proportion, max_m)
            d_b_s.append(d_b)

        draw(
            range(1, parameters.NUMBER_OF_FRAGILE_EDGES),
            d_b_s,
            "x",
            "d/b",
            "d/b depends on x for" + s,
            "fragile_breakage_model/plots/d_b/" + s,
        )


def compare_empirical_and_analytical():
    parameter = parameters.PROBABILITIES_WITH_ALPHA[6]
    string_parameters, p_aa, p_bb, a_type_edges_proportion = parameter
    find_true_evolution_dist(string_parameters, p_aa, p_bb, a_type_edges_proportion, 30)


if __name__ == "__main__":
    compare_empirical_and_analytical()
    # check_monotone_d_divide_b()
    # check_monotone_d_divide_b()
