import time
import os

import numpy as np

import parameters
from aggregate_cycles_info import sum_cycles_info
from utils import (
    CyclesInfo,
    log_experiments,
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


def is_in_cycle(x, edges, y):
    cur = edges[x].to

    if y == cur:
        return True

    while cur != x:
        if cur % 2 == 0:
            cur = cur + 1
        else:
            cur = cur - 1

        if y == cur:
            return True

        if cur == x:
            break

        cur = edges[cur].to

        if y == cur:
            return True

    return False


def set_color(x, edges, colors, color):
    colors[x] = color
    cur = edges[x].to

    while cur != x:
        colors[cur] = color

        if cur % 2 == 0:
            cur = cur + 1
        else:
            cur = cur - 1

        if cur == x:
            break

        colors[cur] = color
        cur = edges[cur].to


def update_cycles(edge1, edge2, joining_type, edges, colors, last_color):
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

    if colors[x] == colors[u]:
        if is_in_cycle(x, edges, y):
            new_color = last_color
        else:
            set_color(x, edges, colors, last_color)
            new_color = last_color + 1
    else:
        set_color(x, edges, colors, colors[x])
        new_color = last_color

    return new_color, new_edges


def compute_cycles_info(colors, edges, max_cycle_len_with_types):
    cycles = {}
    for i, c in enumerate(colors):
        if c in cycles.keys():
            cycles[c].append(edges[i].edge_type)
        else:
            cycles[c] = [edges[i].edge_type]

    cycle_types = {}
    cycles_m = {}
    a_in_non_trivial_cycles = 0
    b_in_non_trivial_cycles = 0
    a_in_non_trivial_cycles_part = 0
    b_in_non_trivial_cycles_part = 0

    for c in cycles.keys():
        # Since we count each edge twice, we should divide it by 2.
        cur_len = len(cycles[c]) // 2
        if cur_len < max_cycle_len_with_types:
            cycle_type = to_cycle_type(cycles[c])
            if cycle_type in cycle_types:
                cycle_types[cycle_type] += 1
            else:
                cycle_types[cycle_type] = 1

        if str(cur_len) in cycles_m:
            cycles_m[str(cur_len)] += 1
        else:
            cycles_m[str(cur_len)] = 1

        if cur_len > 1:
            cnt_a, cnt_b = count_a_b(cycles[c])
            a_in_non_trivial_cycles += cnt_a
            b_in_non_trivial_cycles += cnt_b

            if cur_len < parameters.PART:
                a_in_non_trivial_cycles_part += cnt_a
                b_in_non_trivial_cycles_part += cnt_b

    return CyclesInfo(
        cycle_types,
        cycles_m,
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
    )


def to_cycle_type(edges_type):
    edges_type.sort()
    cycle_type = ""
    for i in range(0, len(edges_type), 2):
        cycle_type += edges_type[i]
    return cycle_type


def count_a_b(edges_type):
    cnt_a = 0
    cnt_b = 0
    for t in edges_type:
        if t == "A":
            cnt_a += 1
        else:
            cnt_b += 1
    return cnt_a // 2, cnt_b // 2


def markov_process(
    n, k, p_aa, p_bb, p_ab, a_type, b_type, edges, max_cycle_len_with_types
):
    steps_cycles_info = []
    len_a = len(a_type)
    len_b = len(b_type)

    probabilities = [p_aa, p_bb, p_ab]

    colors = [0] * (2 * n)
    for i in range(n):
        colors[2 * i] = i
        colors[2 * i + 1] = i

    new_color = n

    steps_cycles_info.append(
        compute_cycles_info(colors, edges, max_cycle_len_with_types)
    )

    for i in range(k):
        dcj_type = np.random.choice(["aa", "bb", "ab"], 1, p=probabilities)[0]
        joining_type = np.random.choice(["xu", "xv"], 1)[0]

        if dcj_type == "aa":
            edge_ind1, edge_ind2 = np.random.choice(len(a_type), 2, replace=False)
            edge1 = a_type[edge_ind1]
            edge2 = a_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            a_type[edge_ind1], a_type[edge_ind2] = new_edges
        elif dcj_type == "bb":
            edge_ind1, edge_ind2 = np.random.choice(len(b_type), 2, replace=False)
            edge1 = b_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            b_type[edge_ind1], b_type[edge_ind2] = new_edges
        else:
            edge_ind1 = np.random.choice(len(a_type), 1)[0]
            edge_ind2 = np.random.choice(len(b_type), 1)[0]

            edge1 = a_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color, new_edges = update_cycles(
                edge1, edge2, joining_type, edges, colors, new_color
            )
            a_type[edge_ind1], b_type[edge_ind2] = new_edges

        steps_cycles_info.append(
            compute_cycles_info(colors, edges, max_cycle_len_with_types)
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
    max_cycle_len_with_types = 6
    max_interesting_cycles_len = parameters.MAX_POSSIBLE_CYCLES_LEN

    for parameter in parameters.PROBABILITIES_WITH_ALPHA[9:12]:
        file = parameter["parameters_str"] + ".csv"
        experiments_in_one_bunch = parameter["experiments_in_one_bunch"]
        experiments_log_path = create_new_directory_for_logging_experiments(experiments_in_one_bunch)

        experiments = []
        print("parameters:", parameter["parameters_str"])

        first_log = True
        for i in range(parameter["number_of_experiments"]):
            a_type, b_type, edges = split_fragile_edges(n, parameter["alpha"])
            experiments.append(
                markov_process(
                    n,
                    parameter["number_of_steps"],
                    parameter["p_aa"],
                    parameter["p_bb"],
                    1 - parameter["p_aa"] - parameter["p_bb"],
                    a_type,
                    b_type,
                    edges,
                    max_cycle_len_with_types,
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
                        max_cycle_len_with_types,
                        max_interesting_cycles_len,
                    ),
                    experiments_log_path + file,
                    open_mode="a",
                    max_cycle_len_with_types=max_cycle_len_with_types,
                    max_interesting_cycles_len=max_interesting_cycles_len,
                )
                experiments = []

    print((time.time() - start_time) / 60 / 60, ".h")


if __name__ == "__main__":
    main()
