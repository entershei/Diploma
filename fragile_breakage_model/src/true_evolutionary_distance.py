import math
import random

import parameters
from generate_directories_names import get_cycles_info_dir
from utils import read_experiments_cycles_info


def compute_alpha(graph):
    a_edges = graph.a_in_non_trivial_cycles + graph.cycle_types["A"]
    b_edges = graph.b_in_non_trivial_cycles + graph.cycle_types["B"]
    error = random.randint(-10, 10)
    a_edges += error
    b_edges += error

    division = a_edges / b_edges
    # division = graph.a_in_non_trivial_cycles / graph.b_in_non_trivial_cycles
    alpha = division / (1 + division)
    return alpha


def compute_probabilities(graph):
    aa_cycles = graph.cycle_types["AA"]
    ab_cycles = graph.cycle_types["AB"]
    bb_cycles = graph.cycle_types["BB"]
    summ = aa_cycles + ab_cycles + bb_cycles
    return aa_cycles / summ, bb_cycles / summ


def compute_probabilities2(graph):
    aaa_cycles = graph.cycle_types["AAA"]
    bbb_cycles = graph.cycle_types["BBB"]
    cycles3 = graph.cycles_m["3"]

    p_aa = math.sqrt(aaa_cycles / cycles3 / 3)
    p_bb = math.sqrt(bbb_cycles / cycles3 / 3)

    p_ab = 1 - p_aa - p_bb

    check = abs(p_ab ** 2 + 2 * p_ab * p_aa - graph.cycle_types["AAB"] / cycles3)

    return p_aa, p_bb


def main():
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = 30

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[5:]:
        file, p_aa, p_bb, alpha = parameter

        print(file)

        graphs = read_experiments_cycles_info(
            get_cycles_info_dir() + file + ".csv",
            max_cycle_len_with_types,
            max_interesting_cycles_len,
            False,
        )[0]

        # Go through each step.
        for k, graph in enumerate(graphs):
            if k < 1:
                continue
            # computed_p_aa, computed_p_bb = compute_probabilities(graph)
            # computed_p_aa2, computed_p_bb2 = compute_probabilities2(graph)

            computed_alpha = compute_alpha(graph)
            print("k:", k)
            print(computed_alpha, abs(alpha - computed_alpha))
            # print(computed_p_aa, computed_p_bb, abs(computed_p_aa - p_aa), abs(computed_p_bb - p_bb))
            # print(computed_p_aa2, computed_p_bb2, abs(computed_p_aa2 - p_aa), abs(computed_p_bb2 - p_bb))

            # assert abs(computed_p_aa - p_aa) < eps
            # assert abs(computed_p_bb - p_bb) < eps


if __name__ == "__main__":
    main()
