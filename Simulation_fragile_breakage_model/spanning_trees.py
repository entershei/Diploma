import csv
import itertools
import time
from functools import cmp_to_key

import numpy as np


def generate_all_dandelion_codes(code, code_len, vertices):
    if code_len < 0:
        return []

    if len(code) == code_len:
        return [code]

    res = []
    for i in vertices:
        res += generate_all_dandelion_codes(code + [i], code_len, vertices)
    return res


def get_cycle(v, connected_to, color, cycle, n):
    if v == 0 or v == n - 1:
        return []

    color[v] = 1
    cycle.append(v)

    if color[connected_to[v]] == 1:
        cycle_beginning = cycle.index(connected_to[v])
        color[v] = 2
        return cycle[cycle_beginning:]
    elif color[connected_to[v]] == 0:
        cycle = get_cycle(connected_to[v], connected_to, color, cycle, n)
        color[v] = 2
        return cycle
    else:
        color[v] = 2
        return []


def rotate(cycle):
    ind = np.argmin(cycle) + 1
    return cycle[ind:] + cycle[:ind]


def find_cycles(code, n):
    cycles = []
    color = [0] * n
    for i in range(1, n - 1):
        if color[i] == 0:
            cycle = get_cycle(i, code, color, [], n)
            if len(cycle) > 0:
                cycles.append(rotate(cycle))
    cycles.sort(key=cmp_to_key(lambda cycle1, cycle2: cycle1[-1] - cycle2[-1]))
    return cycles


def decode_dandelion_code(code):
    code = [-1] + code + [-1]
    n = len(code)
    cycles = find_cycles(code, n)

    used_vertices = [False] * n
    used_vertices[0] = True
    chain = [0] + list(itertools.chain.from_iterable(cycles)) + [n - 1]
    tree = []
    for i in range(1, len(chain)):
        tree.append([chain[i - 1], chain[i]])
        used_vertices[chain[i]] = True

    for i in range(1, n - 1):
        if not used_vertices[i]:
            tree.append([i, code[i]])

    return tree


def convert_to_ab_graph(edges, a_vertices, b_vertices):
    vertices = ["A"] * a_vertices + ["B"] * b_vertices
    for i in range(a_vertices + b_vertices):
        vertices[i] = str(i) + "_" + vertices[i]

    tree = []
    for edge in edges:
        tree.append([vertices[edge[0]], vertices[edge[1]]])
    return tree


def generate_all_spanning_trees(a_vertices, b_vertices, is_special_codes, k1, k3):
    if is_special_codes:
        codes = generate_special_codes(a_vertices, b_vertices, k1, k3)
    else:
        codes = generate_all_dandelion_codes(
            [], a_vertices + b_vertices - 2, range(a_vertices + b_vertices)
        )

    spanning_trees = []
    for code in codes:
        spanning_trees.append(
            convert_to_ab_graph(decode_dandelion_code(code), a_vertices, b_vertices)
        )

    return spanning_trees


def count_edge_types(tree):
    n_aa, n_ab, n_bb = 0, 0, 0
    for edge in tree:
        if edge[0][-1] == "A" and edge[1][-1] == "A":
            n_aa += 1
        elif edge[0][-1] == "B" and edge[1][-1] == "B":
            n_bb += 1
        else:
            n_ab += 1

    return n_aa, n_ab, n_bb


def count_different_trees(spanning_trees):
    different_trees = {}
    for tree in spanning_trees:
        n_aa, n_ab, n_bb = count_edge_types(tree)
        if (n_aa, n_ab, n_bb) in different_trees:
            different_trees[(n_aa, n_ab, n_bb)] += 1
        else:
            different_trees[(n_aa, n_ab, n_bb)] = 1
    return different_trees


def cmp(n1, n2):
    n1_aa, n1_ab, n1_bb = n1
    n2_aa, n2_ab, n2_bb = n2
    if n1_aa != n2_aa:
        return n1_aa - n2_aa
    elif n1_ab != n2_ab:
        return n1_ab - n2_ab
    return n1_bb < n2_bb


def write_log(different_trees, f_name):
    sorted_tree_types = sorted(list(different_trees.keys()), key=cmp_to_key(cmp))
    f_name += ".csv"
    with open(f_name, "w", newline="") as log:
        fieldnames = ["n_aa", "n_ab", "n_bb", "cnt"]
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for trees in sorted_tree_types:
            n_aa, n_ab, n_bb = trees
            log_lens.writerow(
                {
                    "n_aa": n_aa,
                    "n_ab": n_ab,
                    "n_bb": n_bb,
                    "cnt": different_trees[trees],
                }
            )


def log_spanning_trees_types(
        a_vertices, b_vertices, is_special_codes=False, k1=0, k3=0
):
    spanning_trees = generate_all_spanning_trees(
        a_vertices, b_vertices, is_special_codes, k1, k3
    )

    different_trees = count_different_trees(spanning_trees)

    sorted_tree_types = sorted(list(different_trees.keys()), key=cmp_to_key(cmp))
    sum_trees = 0
    for trees in sorted_tree_types:
        sum_trees += different_trees[trees]

    if a_vertices + b_vertices < 2:
        assert len(spanning_trees) == sum_trees == 0
    else:
        if not is_special_codes:
            assert (
                    len(spanning_trees)
                    == sum_trees
                    == (a_vertices + b_vertices) ** (a_vertices + b_vertices - 2)
            )
        else:
            assert len(different_trees) == 1
            assert (
                    len(spanning_trees)
                    == sum_trees
                    == a_vertices ** (k1 + b_vertices - k3 - 1) * b_vertices ** (k3 + a_vertices - k1 - 1)
            )
            assert (sorted_tree_types[0] == (k1, a_vertices + b_vertices - 1 - k1 - k3, k3))

    f_folder = "logs/spanning_trees/"
    f_name = "a_vertices=" + str(a_vertices) + ", b_vertices=" + str(b_vertices)

    if is_special_codes:
        f_folder += "special_trees/"
        f_name += ", k1=" + str(k1) + ", k3=" + str(k3)
    else:
        f_folder += "complete_graph/"

    write_log(different_trees, f_folder + f_name)


def log_different_spanning_trees():
    n = 5
    for a_vertices in range(n):
        for b_vertices in range(n):
            log_spanning_trees_types(a_vertices, b_vertices)


# [2..k1+1]       : [0; n)
# [k1+2..n]       : [n; n+m)
# [n+1..n+k3+1]   : [n, n+m)
# [n+k3+2..n+m-1] : [0; n)
def generate_special_codes(n, m, k1, k3):
    k2 = n - k1 - 1
    k4 = m - k3 - 1
    assert k1 >= 0 and k2 >= 0 and k3 >= 0 and k4 >= 0
    part1 = generate_all_dandelion_codes([], k1, range(n))
    part2 = generate_all_dandelion_codes([], k2, range(n, n + m))
    part3 = generate_all_dandelion_codes([], k3, range(n, n + m))
    part4 = generate_all_dandelion_codes([], k4, range(n))

    codes = []
    for p1 in part1:
        for p2 in part2:
            for p3 in part3:
                for p4 in part4:
                    codes.append(p1 + p2 + p3 + p4)
    return codes


def log_special_spanning_trees():
    n = 5
    for a_vertices in range(1, n):
        for b_vertices in range(1, n):
            for k1 in range(a_vertices):
                for k3 in range(b_vertices):
                    log_spanning_trees_types(a_vertices, b_vertices, True, k1, k3)


def main():
    start_time = time.time()
    # log_different_spanning_trees()
    log_special_spanning_trees()
    print(time.time() - start_time)


if __name__ == "__main__":
    main()
