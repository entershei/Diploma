from results_comparison import compute_analytical_cycles_m

import parameters
from draw_plots import draw, build_parameters_for_plot_title


# Normalized min distance
def compute_analytically_d_n(x, p_aa, p_bb, alpha, max_m):
    dist = 1.0

    for m in range(1, max_m):
        dist -= compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)["all"]

    return dist


# For non trivial cycles normalized sum of cycles lens
def compute_analytically_b_n(x, p_aa, p_bb, alpha):
    return 1 - compute_analytical_cycles_m(1, x, p_aa, p_bb, alpha)["all"]


def compute_d_divide_b(x, p_aa, p_bb, alpha, max_m):
    analytical_min_d = compute_analytically_d_n(x, p_aa, p_bb, alpha, max_m)
    b = compute_analytically_b_n(x, p_aa, p_bb, alpha)
    if b > 0:
        return analytical_min_d / b
    return 0


def draw_d_divide_b():
    for parameter in parameters.PROBABILITIES_WITH_ALPHA:
        s, p_aa, p_bb, alpha = parameter
        print(s)

        max_m = 85
        d_b_s = []
        xs = []
        for k in range(1, parameters.NUMBER_OF_FRAGILE_EDGES):
            x = k / parameters.NUMBER_OF_FRAGILE_EDGES
            xs.append(x)
            d_b = compute_d_divide_b(x, p_aa, p_bb, alpha, max_m)
            d_b_s.append(d_b)

        draw(
            xs,
            d_b_s,
            "x",
            "d/b",
            "d/b depends on x\n" + build_parameters_for_plot_title(p_aa, p_bb, alpha),
            "3d_fragile_breakage_model/plots/d_b/" + s,
        )


if __name__ == "__main__":
    draw_d_divide_b()
