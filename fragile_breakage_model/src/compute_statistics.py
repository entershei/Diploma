from results_comparison import compute_analytical_cycles_m
from generate_directories_names import get_cycles_info_dir

import parameters
from aggregate_cycles_info import read_experiments_cycles_info
from draw_plots import draw


def compute_analytically_normalized_min_dist(x, p_aa, p_bb, alpha, max_m):
    dist = 1.0

    for m in range(1, max_m):
        dist -= compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]

    return dist


def compute_analytically_normalized_number_non_trivial_cycles(
    x, p_aa, p_bb, alpha, max_m
):
    # b = 0.0
    # for m in range(2, max_m):
    #     b += m * compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]
    # return b
    return 1 - compute_analytical_cycles_m(1, x, p_aa, p_bb, alpha)["all"]


def compute_d_with_different_parameters(max_m):
    p = list(map(lambda a: a / 10, range(1, 10)))
    for p_aa in p:
        for p_bb in p:
            if p_aa + p_bb < 1:
                for x in p:
                    for alpha in p:
                        d, iterations = compute_analytically_normalized_min_dist(
                            x, p_aa, p_bb, alpha, max_m=max_m
                        )
                        print(d, iterations, p_aa, p_bb, x, alpha)

    # check_d_b(get_parameters_as_string() + "paa0,5_pbb0,45_alpha0,5_10000_experiments.csv")


def compare_min_d(parameters_str, p_aa, p_bb, a_type_edges_proportion, max_m):
    aggregated_info = read_experiments_cycles_info(
        get_cycles_info_dir() + parameters_str + ".csv",
        5,
        parameters.MAX_POSSIBLE_CYCLES_LEN,
        False,
    )[0]

    for k, info in enumerate(aggregated_info):
        empirical_min_d = 0
        for cycle_len in info.cycles_m.keys():
            if int(cycle_len) > 1:
                empirical_min_d += info.cycles_m[cycle_len] * (int(cycle_len) - 1)

        empirical_min_d /= parameters.NUMBER_OF_FRAGILE_EDGES
        analytical_min_d = compute_analytically_normalized_min_dist(
            k / parameters.NUMBER_OF_FRAGILE_EDGES,
            p_aa,
            p_bb,
            a_type_edges_proportion,
            max_m,
        )
        # print(
        #     k, abs(empirical_min_d - analytical_min_d), empirical_min_d, analytical_min_d
        # )
        if abs(empirical_min_d - analytical_min_d) <= 1e-2:
            print(
                k,
                abs(empirical_min_d - analytical_min_d),
                empirical_min_d,
                analytical_min_d,
            )


def compute_d_divide_b(x, p_aa, p_bb, alpha, max_m):
    analytical_min_d = compute_analytically_normalized_min_dist(
        x, p_aa, p_bb, alpha, max_m
    )
    b = compute_analytically_normalized_number_non_trivial_cycles(
        x, p_aa, p_bb, alpha, max_m
    )
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


def compare_different_min_d():
    parameter = parameters.PROBABILITIES_WITH_ALPHA[1]
    string_parameters, p_aa, p_bb, a_type_edges_proportion = parameter
    compare_min_d(string_parameters, p_aa, p_bb, a_type_edges_proportion, 111)


if __name__ == "__main__":
    # compute_d_with_different_parameters()
    # compare_different_min_d()
    check_monotone_d_divide_b()
