import csv
import time

import numpy as np


def read_number_of_fragile_regions_and_steps():
    print("Enter number of fragile regions")
    n = int(input())

    print("Enter number of steps")
    k = int(input())

    return n, k


def read_probabilities():
    print("Enter probability of rearrangement between A-type fragile regions:")
    p_aa = float(input())

    print("Enter probability of rearrangement between B-type fragile regions:")
    p_bb = float(input())

    # Probability of rearrangement between A-type and D-type fragile regions.
    p_ab = 1 - p_aa - p_bb

    assert p_aa >= 0 and p_bb >= 0 and p_ab >= 0

    return p_aa, p_bb, p_ab


class EdgeInfo:
    to = 0
    edge_type = "A"

    def __init__(self, to, edge_type):
        self.to = to
        self.edge_type = edge_type


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


class CyclesInfo:
    num_all_cycles = 0
    # cycles with length from 1 to 5 with different types
    num_n_cycles = {}
    max_len = 0

    def __init__(self, num_all_cycles, num_n_cycles, max_len):
        self.num_all_cycles = num_all_cycles
        self.num_n_cycles = num_n_cycles
        self.max_len = max_len


# Возвращает количество циклов; [число циклов длины 1: (с хрупким ребром A; с хрупким ребром B),
# длины 2: (с хрупкими ребрами A, A; с хрупкими ребрами A, B; с хрупкими ребрами B, B),
# длины 3: (с хрупкими ребрами A, A, A; с хрупкими ребрами A, A, B; с хрупкими ребрами A, B, B;
#           с хрупкими ребрами B, B, B),
# длины 4: (с хрупкими ребрами A, A, A, A; с хрупкими ребрами A, A, A, B; с хрупкими ребрами A, A, B, B;
#           с хрупкими ребрами A, B, B, B; с хрупкими ребрами B, B, B, B),
# длины 5: (с хрупкими ребрами A, A, A, A, A; с хрупкими ребрами A, A, A, A, B; с хрупкими ребрами A, A, A, B, B;
#           с хрупкими ребрами A, A, B, B, B; с хрупкими ребрами A, B, B, B, B; с хрупкими ребрами B, B, B, B, B)];
# длина максимального цикла.
# Длина цикла измеряется в количестве хрупких ребер в нем.
def compute_cycles_info(colors, edges):
    cycles = {}
    for i, c in enumerate(colors):
        if c in cycles.keys():
            cycles[c].append(edges[i].edge_type)
        else:
            cycles[c] = [edges[i].edge_type]

    cycles_info = {}
    max_len = 0
    for c in cycles.keys():
        cur_len = len(cycles[c]) // 2
        if cur_len < 6:
            cycle_type = to_cycle_type(cycles[c])
            if cycle_type in cycles_info:
                cycles_info[cycle_type] += 1
            else:
                cycles_info[cycle_type] = 1
        max_len = max(max_len, cur_len)

    return CyclesInfo(len(cycles), cycles_info, max_len)


def to_cycle_type(edges_type):
    edges_type.sort()
    cycle_type = ""
    for i in range(0, len(edges_type), 2):
        cycle_type += edges_type[i]
    return cycle_type


def markov_process(n, k, p_aa, p_bb, p_ab, a_type, b_type, edges):
    steps_cycles_info = []
    len_a = len(a_type)
    len_b = len(b_type)

    probabilities = [p_aa, p_bb, p_ab]

    colors = [0] * (2 * n)
    for i in range(n):
        colors[2 * i] = i
        colors[2 * i + 1] = i

    new_color = n

    steps_cycles_info.append(compute_cycles_info(colors, edges))

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

        steps_cycles_info.append(compute_cycles_info(colors, edges))

        assert len_a == len(a_type) and len_b == len(b_type)

    return steps_cycles_info


def aggregate_cycles_info(experiments):
    num_experiments = len(experiments)
    num_steps_of_markov_process = len(experiments[0])

    aggregated_cycles_info = []
    for i in range(num_steps_of_markov_process):
        average_all_cycles_num = 0
        average_cnt_cycles = {
            "A": 0,
            "B": 0,
            "AA": 0,
            "AB": 0,
            "BB": 0,
            "AAA": 0,
            "AAB": 0,
            "ABB": 0,
            "BBB": 0,
            "AAAA": 0,
            "AAAB": 0,
            "AABB": 0,
            "ABBB": 0,
            "BBBB": 0,
            "AAAAA": 0,
            "AAAAB": 0,
            "AAABB": 0,
            "AABBB": 0,
            "ABBBB": 0,
            "BBBBB": 0,
        }

        average_max_cycles_len = 0
        for j in range(num_experiments):
            average_all_cycles_num += experiments[j][i].num_all_cycles
            for cycles_types in experiments[j][i].num_n_cycles.keys():
                average_cnt_cycles[cycles_types] += experiments[j][i].num_n_cycles[
                    cycles_types
                ]

            average_max_cycles_len += experiments[j][i].max_len

        average_all_cycles_num /= num_experiments
        for cycles_types in average_cnt_cycles.keys():
            average_cnt_cycles[cycles_types] /= num_experiments
        average_max_cycles_len /= num_experiments

        aggregated_cycles_info.append(
            CyclesInfo(
                average_all_cycles_num, average_cnt_cycles, average_max_cycles_len
            )
        )

    return aggregated_cycles_info


def log_aggregated_results(n, aggregated_cycles_info, f):
    with open(f, "w", newline="") as f_log_lens:
        fieldnames = [
            "n",
            "k",
            "all-cycles",
            "A-cycles",
            "B-cycles",
            "AA-cycles",
            "AB-cycles",
            "BB-cycles",
            "AAA-cycles",
            "AAB-cycles",
            "ABB-cycles",
            "BBB-cycles",
            "AAAA-cycles",
            "AAAB-cycles",
            "AABB-cycles",
            "ABBB-cycles",
            "BBBB-cycles",
            "AAAAA-cycles",
            "AAAAB-cycles",
            "AAABB-cycles",
            "AABBB-cycles",
            "ABBBB-cycles",
            "BBBBB-cycles",
            "max_cycle_len",
        ]
        log_cycles_info = csv.DictWriter(f_log_lens, fieldnames=fieldnames)
        log_cycles_info.writeheader()

        for step, info in enumerate(aggregated_cycles_info):
            log_cycles_info.writerow(
                {
                    "n": n,
                    "k": step,
                    "all-cycles": info.num_all_cycles,
                    "A-cycles": info.num_n_cycles["A"],
                    "B-cycles": info.num_n_cycles["B"],
                    "AA-cycles": info.num_n_cycles["AA"],
                    "AB-cycles": info.num_n_cycles["AB"],
                    "BB-cycles": info.num_n_cycles["BB"],
                    "AAA-cycles": info.num_n_cycles["AAA"],
                    "AAB-cycles": info.num_n_cycles["AAB"],
                    "ABB-cycles": info.num_n_cycles["ABB"],
                    "BBB-cycles": info.num_n_cycles["BBB"],
                    "AAAA-cycles": info.num_n_cycles["AAAA"],
                    "AAAB-cycles": info.num_n_cycles["AAAB"],
                    "AABB-cycles": info.num_n_cycles["AABB"],
                    "ABBB-cycles": info.num_n_cycles["ABBB"],
                    "BBBB-cycles": info.num_n_cycles["BBBB"],
                    "AAAAA-cycles": info.num_n_cycles["AAAAA"],
                    "AAAAB-cycles": info.num_n_cycles["AAAAB"],
                    "AAABB-cycles": info.num_n_cycles["AAABB"],
                    "AABBB-cycles": info.num_n_cycles["AABBB"],
                    "ABBBB-cycles": info.num_n_cycles["ABBBB"],
                    "BBBBB-cycles": info.num_n_cycles["BBBBB"],
                    "max_cycle_len": info.max_len,
                }
            )


def main():
    # n, k = read_number_of_regions_and_steps()
    # p_aa, p_bb, p_ab = read_probabilities()
    start_time = time.time()
    n, p_aa, p_bb = 1000, 0.45, 0.45
    a_type_edges_proportion = 0.5

    p_ab = 1 - p_aa - p_bb
    k = n * 2

    experiments = []
    different_fragile_edges_splits = 50
    markov_process_experiments = 20
    for i in range(different_fragile_edges_splits):
        a_type, b_type, edges = split_fragile_edges(n, a_type_edges_proportion)
        for j in range(markov_process_experiments):
            experiments.append(
                markov_process(
                    n, k, p_aa, p_bb, p_ab, a_type.copy(), b_type.copy(), edges.copy()
                )
            )

    aggregated_cycles_info = aggregate_cycles_info(experiments)

    log_aggregated_results(
        n,
        aggregated_cycles_info,
        "logs/cycles_info/n1000/"
        + str(different_fragile_edges_splits)
        + "_"
        + str(markov_process_experiments)
        + "_experiments/"
        + "paa0_45_pbb0_45_alpha0_5.csv",
    )
    print(time.time() - start_time)


if __name__ == "__main__":
    main()
