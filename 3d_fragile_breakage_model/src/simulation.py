import time
import os

import numpy as np

from collections import defaultdict

import parameters
from aggregate_cycles_info import sum_cycles_info
from compute_statistics import compute_d_by_cycles, compute_b_by_cycles
from utils import (
    CyclesInfo,
    log_experiments,
    define_cycles_representative,
)
from generate_directories_names import create_new_directory_for_logging_experiments


def read_parameters_from_console():
    print("Enter number of fragile regions")
    n = int(input())

    print("Enter number of steps")
    k = int(input())

    print("Enter probability of rearrangement between A-type fragile regions:")
    p_aa = float(input())

    print("Enter probability of rearrangement between B-type fragile regions:")
    p_bb = float(input())

    # Probability of rearrangement between A-type and D-type fragile regions.
    p_ab = 1 - p_aa - p_bb

    assert p_aa >= 0 and p_bb >= 0 and p_ab >= 0

    return n, k, p_aa, p_bb, p_ab


class EdgeInfo:
    to = 0
    edge_type = "A"

    def __init__(self, to, edge_type):
        self.to = to
        self.edge_type = edge_type


# Split only edges for Q genome.
def split_fragile_edges(n, a_type_edges_proportion):
    a_type = []
    b_type = []
    edges = []

    will_be_a_type = np.random.choice(
        range(0, n), int(a_type_edges_proportion * n), replace=False
    )
    is_a_type = [False] * n
    for e in will_be_a_type:
        is_a_type[e] = True

    for i in range(n):
        edge = [2 * i, 2 * i + 1]
        if is_a_type[i]:
            a_type.append(edge)
            edges.append(EdgeInfo(2 * i + 1, "A"))
            edges.append(EdgeInfo(2 * i, "A"))
        else:
            b_type.append(edge)
            edges.append(EdgeInfo(2 * i + 1, "B"))
            edges.append(EdgeInfo(2 * i, "B"))

    return a_type, b_type, edges


def update_cycles(edge1, edge2, joining_type, edges):
    x, y = edge1
    u, v = edge2

    if joining_type == "xu":
        new_edge1_type, new_edge2_type = np.random.choice(
            [edges[x].edge_type, edges[u].edge_type], 2, replace=False
        )

        edges[x], edges[u] = EdgeInfo(u, new_edge1_type), EdgeInfo(x, new_edge1_type)
        edges[y], edges[v] = EdgeInfo(v, new_edge2_type), EdgeInfo(y, new_edge2_type)

        new_edges = [[x, u], [y, v]]
        if new_edge1_type != new_edge2_type and new_edge1_type == "B":
            new_edges.reverse()
    else:
        new_edge1_type, new_edge2_type = np.random.choice(
            [edges[x].edge_type, edges[v].edge_type], 2, replace=False
        )

        edges[x], edges[v] = EdgeInfo(v, new_edge1_type), EdgeInfo(x, new_edge1_type)
        edges[y], edges[u] = EdgeInfo(u, new_edge2_type), EdgeInfo(y, new_edge2_type)

        new_edges = [[x, v], [y, u]]
        if new_edge1_type != new_edge2_type and new_edge1_type == "B":
            new_edges.reverse()

    return new_edges


def get_cycle(start_v, used, edges):
    used[start_v] = True

    cycle_edges_types = [edges[start_v].edge_type]

    cur_v = edges[start_v].to
    while cur_v != start_v:
        used[cur_v] = True

        if cur_v % 2 == 0:
            cur_v = cur_v + 1
        else:
            cur_v = cur_v - 1

        if cur_v == start_v:
            break

        used[cur_v] = True
        cycle_edges_types.append(edges[cur_v].edge_type)
        cur_v = edges[cur_v].to

    return "".join(cycle_edges_types)


def compute_cycles_info(edges, max_cycle_len_with_types, cycle_to_represent):
    cycles_with_edges_order = defaultdict(int)
    cycles_m = defaultdict(int)
    a_in_non_trivial_cycles = 0
    b_in_non_trivial_cycles = 0
    a_in_non_trivial_cycles_part = 0
    b_in_non_trivial_cycles_part = 0

    vertices = len(edges)
    used = [False] * len(edges)
    for v in range(vertices):
        if not used[v]:
            cycle_edges = get_cycle(v, used, edges)
            cycle_len = len(cycle_edges)

            if cycle_len < max_cycle_len_with_types:
                cycles_with_edges_order[cycle_to_represent[cycle_edges]] += 1

            cycles_m[str(cycle_len)] += 1

            if cycle_len > 1:
                a_in_non_trivial_cycles += cycle_edges.count("A")
                b_in_non_trivial_cycles += cycle_edges.count("B")

                if cycle_len < parameters.PART:
                    a_in_non_trivial_cycles_part += cycle_edges.count("A")
                    b_in_non_trivial_cycles_part += cycle_edges.count("B")

    return CyclesInfo(
        dict(cycles_m),
        dict(cycles_with_edges_order),
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
        compute_d_by_cycles(cycles_m),
        compute_b_by_cycles(cycles_m),
    )


def markov_process(
    k,
    p_aa,
    p_bb,
    p_ab,
    a_type,
    b_type,
    edges,
    max_cycle_len_with_types,
    cycle_to_represent,
):
    steps_cycles_info = []
    len_a = len(a_type)
    len_b = len(b_type)

    probabilities = [p_aa, p_bb, p_ab]

    steps_cycles_info.append(
        compute_cycles_info(edges, max_cycle_len_with_types, cycle_to_represent)
    )

    for i in range(k):
        dcj_type = np.random.choice(["aa", "bb", "ab"], 1, p=probabilities)[0]
        joining_type = np.random.choice(["xu", "xv"], 1)[0]

        if dcj_type == "aa":
            edge_ind1, edge_ind2 = np.random.choice(len(a_type), 2, replace=False)
            edge1 = a_type[edge_ind1]
            edge2 = a_type[edge_ind2]

            new_edges = update_cycles(edge1, edge2, joining_type, edges)
            a_type[edge_ind1], a_type[edge_ind2] = new_edges
        elif dcj_type == "bb":
            edge_ind1, edge_ind2 = np.random.choice(len(b_type), 2, replace=False)
            edge1 = b_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_edges = update_cycles(edge1, edge2, joining_type, edges)
            b_type[edge_ind1], b_type[edge_ind2] = new_edges
        else:
            edge_ind1 = np.random.choice(len(a_type), 1)[0]
            edge_ind2 = np.random.choice(len(b_type), 1)[0]

            edge1 = a_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_edges = update_cycles(edge1, edge2, joining_type, edges)
            a_type[edge_ind1], b_type[edge_ind2] = new_edges

        steps_cycles_info.append(
            compute_cycles_info(edges, max_cycle_len_with_types, cycle_to_represent)
        )

        assert len_a == len(a_type) and len_b == len(b_type)

    return steps_cycles_info


def remove_previous_log(file):
    if os.path.exists(file):
        os.remove(file)


def main():
    # n, k, p_aa, p_bb, p_ab = read_parameters_from_console()
    start_time = time.time()

    n = parameters.NUMBER_OF_FRAGILE_EDGES
    max_cycle_len_with_types = 7
    max_interesting_cycles_len = parameters.MAX_POSSIBLE_CYCLES_LEN

    to_represent, representatives = define_cycles_representative(
        max_cycle_len_with_types
    )

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[15:]:
        file = parameter["parameters_str"] + ".csv"
        experiments_in_one_bunch = parameter["experiments_in_one_bunch"]
        experiments_log_path = create_new_directory_for_logging_experiments(
            experiments_in_one_bunch
        )

        experiments = []
        print("parameters:", parameter["parameters_str"])

        first_log = True
        for i in range(parameter["number_of_experiments"]):
            a_type, b_type, edges = split_fragile_edges(n, parameter["alpha"])
            experiments.append(
                markov_process(
                    parameter["number_of_steps"],
                    parameter["p_aa"],
                    parameter["p_bb"],
                    1 - parameter["p_aa"] - parameter["p_bb"],
                    a_type,
                    b_type,
                    edges,
                    max_cycle_len_with_types,
                    to_represent,
                )
            )

            if len(experiments) == experiments_in_one_bunch:
                if first_log:
                    remove_previous_log(experiments_log_path + file)
                    first_log = False

                # Записываем сумму результатов экспериментов
                if i % 100 == 0:
                    print("i:", i, ", time: ", (time.time() - start_time) / 60, " m.")
                log_experiments(
                    sum_cycles_info(
                        experiments,
                        max_interesting_cycles_len,
                        representatives,
                    ),
                    experiments_log_path + file,
                    open_mode="a",
                    max_interesting_cycles_len=max_interesting_cycles_len,
                    cycles_representatives=representatives,
                )
                experiments = []

    print((time.time() - start_time) / 60 / 60, ".h")


if __name__ == "__main__":
    main()
