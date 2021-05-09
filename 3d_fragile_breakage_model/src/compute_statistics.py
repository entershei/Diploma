from results_comparison import compute_analytical_cycles_m

import parameters
from draw_plots import draw, build_parameters_for_plot_title, draw_plots

from src.generate_directories_names import get_cycles_info_dir
from src.utils import read_experiments_cycles_info
from pathlib import Path


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


def for_finding_parameters():
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = 6

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[:4]:
        string_parameters, p_aa, p_bb, alpha = parameter
        file = get_cycles_info_dir() + string_parameters + ".csv"
        p_ab = 1 - p_aa - p_bb
        beta = 1 - alpha

        Path(
            "3d_fragile_breakage_model/plots/statistic/" + string_parameters + "/"
        ).mkdir(parents=True, exist_ok=True)

        cycles_info = read_experiments_cycles_info(
            file, max_cycle_len_with_types, max_interesting_cycles_len, False
        )[0][1:]

        c_aa_divide_c2 = list(
            map(
                lambda cycle_info: cycle_info.cycle_types["AA"]
                / cycle_info.cycles_m["2"],
                cycles_info,
            )
        )
        c_ab_divide_c2 = list(
            map(
                lambda cycle_info: cycle_info.cycle_types["AB"]
                / cycle_info.cycles_m["2"],
                cycles_info,
            )
        )
        c_bb_divide_c2 = list(
            map(
                lambda cycle_info: cycle_info.cycle_types["BB"]
                / cycle_info.cycles_m["2"],
                cycles_info,
            )
        )
        c_aa_divide_c_bb = list(
            map(
                lambda cycle_info: cycle_info.cycle_types["AA"]
                / cycle_info.cycle_types["BB"],
                cycles_info,
            )
        )

        a_edges_divide_non_trivial = list(
            map(
                lambda cycle_info: cycle_info.a_in_non_trivial_cycles
                / (
                    cycle_info.a_in_non_trivial_cycles
                    + cycle_info.b_in_non_trivial_cycles
                ),
                cycles_info,
            )
        )
        b_edges_divide_non_trivial = list(
            map(
                lambda cycle_info: cycle_info.b_in_non_trivial_cycles
                / (
                    cycle_info.a_in_non_trivial_cycles
                    + cycle_info.b_in_non_trivial_cycles
                ),
                cycles_info,
            )
        )

        a_edges_divide_non_trivial_part = list(
            map(
                lambda cycle_info: cycle_info.a_in_non_trivial_cycles_part
                / (
                    cycle_info.a_in_non_trivial_cycles_part
                    + cycle_info.b_in_non_trivial_cycles_part
                ),
                cycles_info,
            )
        )
        b_edges_divide_non_trivial_part = list(
            map(
                lambda cycle_info: cycle_info.b_in_non_trivial_cycles_part
                / (
                    cycle_info.a_in_non_trivial_cycles_part
                    + cycle_info.b_in_non_trivial_cycles_part
                ),
                cycles_info,
            )
        )

        parameters_for_plot_title = build_parameters_for_plot_title(p_aa, p_bb, alpha)
        steps = len(cycles_info)
        xs = range(1, steps + 1)
        save_path = (
            "3d_fragile_breakage_model/plots/statistic/" + string_parameters + "/"
        )

        draw_plots(
            xs,
            [
                {"plot": c_aa_divide_c2, "label": "c_aa/c2", "color": "red"},
                {"plot": [p_aa] * steps, "label": "p_aa", "color": "black"},
                {"plot": [alpha] * steps, "label": "alpha", "color": "blue"},
            ],
            "steps",
            "normalized cycles",
            "c_aa/c2 depends on steps,\n" + parameters_for_plot_title,
            save_path + "c_aa_c2",
        )
        draw_plots(
            xs,
            [
                {"plot": c_ab_divide_c2, "label": "c_ab/c2", "color": "red"},
                {"plot": [p_ab] * steps, "label": "p_ab", "color": "black"},
                {"plot": [alpha] * steps, "label": "alpha", "color": "blue"},
                {"plot": [beta] * steps, "label": "beta", "color": "slateblue"},
            ],
            "steps",
            "normalized cycles",
            "c_ab/c2 depends on steps,\n" + parameters_for_plot_title,
            save_path + "c_ab_c2",
        )
        draw_plots(
            xs,
            [
                {"plot": c_bb_divide_c2, "label": "c_bb/c2", "color": "red"},
                {"plot": [p_bb] * steps, "label": "p_bb", "color": "black"},
                {"plot": [alpha] * steps, "label": "alpha", "color": "blue"},
            ],
            "steps",
            "normalized cycles",
            "c_bb/c2 depends on steps,\n" + parameters_for_plot_title,
            save_path + "c_bb_c2",
        )
        draw_plots(
            xs,
            [
                {"plot": c_aa_divide_c_bb, "label": "c_aa/c_bb", "color": "red"},
                {"plot": [p_aa / p_bb] * steps, "label": "p_aa", "color": "black"},
            ],
            "steps",
            "normalized cycles",
            "c_aa/c_bb depends on steps,\n" + parameters_for_plot_title,
            save_path + "c_aa_c_bb",
        )

        draw_plots(
            xs,
            [
                {
                    "plot": a_edges_divide_non_trivial,
                    "label": "A/(A+B) edges",
                    "color": "red",
                },
                {"plot": [alpha] * steps, "label": "alpha", "color": "blue"},
            ],
            "steps",
            "normalized edges",
            "A-edges/(A-edges + B-edges) in non trivial cycles,\n"
            + parameters_for_plot_title,
            save_path + "a_edges",
        )
        draw_plots(
            xs,
            [
                {
                    "plot": b_edges_divide_non_trivial,
                    "label": "B/(A+B) edges",
                    "color": "red",
                },
                {"plot": [beta] * steps, "label": "beta", "color": "blue"},
            ],
            "steps",
            "normalized edges",
            "B-edges/(A-edges + B-edges) in non trivial cycles,\n"
            + parameters_for_plot_title,
            save_path + "b_edges",
        )

        draw_plots(
            xs,
            [
                {
                    "plot": a_edges_divide_non_trivial_part,
                    "label": "A/(A+B) part edges",
                    "color": "red",
                },
                {"plot": [alpha] * steps, "label": "alpha", "color": "blue"},
            ],
            "steps",
            "normalized edges",
            "A-edges/(A-edges + B-edges) in non trivial cycles with len <= 50,\n"
            + parameters_for_plot_title,
            save_path + "a_edges_part",
        )
        draw_plots(
            xs,
            [
                {
                    "plot": b_edges_divide_non_trivial_part,
                    "label": "B/(A+B) part edges",
                    "color": "red",
                },
                {"plot": [beta] * steps, "label": "beta", "color": "blue"},
            ],
            "steps",
            "normalized edges",
            "B-edges/(A-edges + B-edges) in non trivial cycles with len <= 50,\n"
            + parameters_for_plot_title,
            save_path + "b_edges_part",
        )


if __name__ == "__main__":
    # draw_d_divide_b()
    for_finding_parameters()
