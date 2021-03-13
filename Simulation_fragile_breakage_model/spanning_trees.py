import csv
import itertools
from functools import cmp_to_key

import numpy as np


def generate_all_dandelion_codes(code, n):
    if n < 2:
        return []

    if len(code) == n - 2:
        return [code]

    res = []
    for i in range(n):
        res += generate_all_dandelion_codes(code + [i], n)
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


def convert_to_ab_graph(edges, vertices):
    tree = []
    for edge in edges:
        tree.append([vertices[edge[0]], vertices[edge[1]]])
    return tree


def generate_all_spanning_trees(a_vertices, b_vertices):
    codes = generate_all_dandelion_codes([], a_vertices + b_vertices)

    vertices = ["A"] * a_vertices + ["B"] * b_vertices
    for i in range(a_vertices + b_vertices):
        vertices[i] = str(i) + "_" + vertices[i]

    spanning_trees = []
    for code in codes:
        spanning_trees.append(
            convert_to_ab_graph(decode_dandelion_code(code), vertices)
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


def log_spanning_trees_types(a_vertices, b_vertices):
    spanning_trees = generate_all_spanning_trees(a_vertices, b_vertices)
    different_trees = count_different_trees(spanning_trees)

    sorted_tree_types = sorted(list(different_trees.keys()), key=cmp_to_key(cmp))
    sum_trees = 0
    for trees in sorted_tree_types:
        sum_trees += different_trees[trees]

    if a_vertices + b_vertices < 2:
        assert sum_trees == 0
    else:
        assert sum_trees == (a_vertices + b_vertices) ** (a_vertices + b_vertices - 2)

    f_folder = "logs/spanning_trees/"
    f_name = "a_vertices=" + str(a_vertices) + ", b_vertices=" + str(b_vertices)
    write_log(different_trees, f_folder + f_name)


def main():
    n = 5
    for a_vertices in range(n):
        for b_vertices in range(n):
            log_spanning_trees_types(a_vertices, b_vertices)


if __name__ == "__main__":
    main()
