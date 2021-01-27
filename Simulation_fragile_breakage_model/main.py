import csv
import random

import matplotlib.pyplot as plt
import numpy as np


def read_number_of_regions_and_steps():
    print('Enter number of fragile regions')
    n = int(input())

    print('Enter number of steps')
    k = int(input())

    return n, k


def read_probabilities():
    print('Enter probability of rearrangement between A-type fragile regions:')
    p_aa = float(input())

    print('Enter probability of rearrangement between B-type fragile regions:')
    p_bb = float(input())

    # Probability of rearrangement between A-type and D-type fragile regions.
    p_ab = 1 - p_aa - p_bb

    assert p_aa >= 0 and p_bb >= 0 and p_ab >= 0

    return p_aa, p_bb, p_ab


def is_in_cycle(x, edges, y):
    cur = edges[x]

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

        cur = edges[cur]

        if y == cur:
            return True

    return False


def set_color(x, edges, colors, color):
    colors[x] = color
    cur = edges[x]

    while cur != x:
        colors[cur] = color

        if cur % 2 == 0:
            cur = cur + 1
        else:
            cur = cur - 1

        if cur == x:
            break

        colors[cur] = color
        cur = edges[cur]


def update_cycles_number(edge1, edge2, joining_type, edges, cnt_cycles, colors, new_color):
    x, y = edge1
    u, v = edge2

    if joining_type == 'xu':
        edges[x], edges[u] = u, x
        edges[y], edges[v] = v, y
    else:
        edges[x], edges[v] = v, x
        edges[y], edges[u] = u, y

    if colors[x] == colors[u]:
        if is_in_cycle(x, edges, y):
            return cnt_cycles, new_color
        else:
            set_color(x, edges, colors, new_color)
            return cnt_cycles + 1, new_color + 1
    else:
        set_color(x, edges, colors, colors[x])
        return cnt_cycles - 1, new_color


def update_edges(edge1, edge2, join_type):
    x, y = edge1
    u, v = edge2

    if join_type == 'xu':
        new_edges = [[x, u], [y, v]]
    else:
        new_edges = [[x, v], [y, u]]

    random.shuffle(new_edges)
    return new_edges


def compute_cycles_lens(n, k, colors, log_lens):
    cycles = {}
    for c in colors:
        if c in cycles.keys():
            cycles[c] += 1
        else:
            cycles[c] = 1

    lens = [0] * 6
    max_len = 0
    for c in cycles.keys():
        cur_len = cycles[c] // 2
        if cur_len < 6:
            lens[cur_len] += 1
        max_len = max(max_len, cur_len)

    log_lens.writerow(
        {'n': n, 'k': k, 'cycles': len(cycles), '1-cycles': lens[1], '2-cycles': lens[2], '3-cycles': lens[3],
         '4-cycles': lens[4], '5-cycles': lens[5], 'max_cycle_len': max_len})


def markov_process(n, k, p_aa, p_bb, p_ab, a_type, b_type, edges, f_log, log_lens):
    len_a = len(a_type)
    len_b = len(b_type)

    probabilities = [p_aa, p_bb, p_ab]

    cycles = [n]
    colors = [0] * (2 * n)
    for i in range(n):
        colors[2 * i] = i
        colors[2 * i + 1] = i

    cnt_cycles = n
    new_color = n

    cnt_aa = 0
    cnt_bb = 0
    cnt_ab = 0

    compute_cycles_lens(n, 0, colors, log_lens)

    for i in range(k):
        dcj_type = np.random.choice(['aa', 'bb', 'ab'], 1, p=probabilities)[0]
        joining_type = np.random.choice(['xu', 'xv'], 1)[0]

        if dcj_type == 'aa':
            cnt_aa += 1
            edge_ind1, edge_ind2 = np.random.choice(len(a_type), 2, replace=False)
            edge1 = a_type[edge_ind1]
            edge2 = a_type[edge_ind2]

            cnt_cycles, new_color = update_cycles_number(edge1, edge2, joining_type, edges, cnt_cycles, colors,
                                                         new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            a_type[edge_ind1], a_type[edge_ind2] = new_edges
        elif dcj_type == 'bb':
            cnt_bb += 1
            edge_ind1, edge_ind2 = np.random.choice(len(b_type), 2, replace=False)
            edge1 = b_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            cnt_cycles, new_color = update_cycles_number(edge1, edge2, joining_type, edges, cnt_cycles, colors,
                                                         new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            b_type[edge_ind1], b_type[edge_ind2] = new_edges
        else:
            cnt_ab += 1
            edge_ind1 = np.random.choice(len(a_type), 1)[0]
            edge_ind2 = np.random.choice(len(b_type), 1)[0]

            edge1 = a_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            cnt_cycles, new_color = update_cycles_number(edge1, edge2, joining_type, edges, cnt_cycles, colors,
                                                         new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            a_type[edge_ind1], b_type[edge_ind2] = new_edges

        cycles.append(cnt_cycles)

        assert len_a == len(a_type) and len_b == len(b_type)

        if i < 10:
            f_log.write(
                'Step ' + str(i) + ':\nRegions = ' + dcj_type + ' ' + joining_type + ' (' + str(edge1[0]) + ', ' + str(
                    edge1[1]) + ') ('
                + str(edge2[0]) + ', ' + str(edge2[1]) + ')    ')
            f_log.write('Number of cycles = ' + str(cnt_cycles) + '\n')

            # f_log.write('Edges:\n')
            # for e in edges:
            #     f_log.write(str(e) + ' ')
            #
            # f_log.write('\nA-type:\n')
            # for e in a_type:
            #     f_log.write('(' + str(e[0]) + ', ' + str(e[1]) + ')  ')
            #
            # f_log.write('\nB-type:\n')
            # for e in b_type:
            #     f_log.write('(' + str(e[0]) + ', ' + str(e[1]) + ')  ')
            #
            # f_log.write('\nColors:\n')
            # for c in colors:
            #     f_log.write(str(c) + ' ')

            f_log.write('\n\n')

        compute_cycles_lens(n, i + 1, colors, log_lens)

    f_log.write('\nNumber aa, bb, ab = ' + str(cnt_aa) + ' ' + str(cnt_bb) + ' ' + str(cnt_ab) + '\n')

    return cycles


def draw_number_of_cycles(cycles):
    k = len(cycles)
    x = list(range(k))
    plt.plot(x, cycles)
    plt.title('Number of cycles depends of number of swaps')
    plt.xlabel('Number of swaps')
    plt.ylabel('Cycles')
    plt.show()


def split_fragile_edges(n, p):
    a_type = []
    b_type = []
    edges = [0] * (2 * n)

    for i in range(n):
        a_or_b = np.random.choice(['a', 'b'], 1, p=p)[0]
        edge = [2 * i, 2 * i + 1]
        if a_or_b == 'a':
            a_type.append(edge)
        else:
            b_type.append(edge)
        edges[2 * i] = 2 * i + 1
        edges[2 * i + 1] = 2 * i

    return a_type, b_type, edges


def main():
    # n, k = read_number_of_regions_and_steps()
    # p_aa, p_bb, p_ab = read_probabilities()
    n, p_aa, p_bb, p_ab = 2000, 0.5, 0.45, 0.05
    k = n * 2

    p_split_a_b = [0.5, 0.5]
    a_type, b_type, edges = split_fragile_edges(n, p_split_a_b)

    f_log = open('log.txt', 'w')

    with open('log_lens.csv', 'w', newline='') as f_log_lens:
        fieldnames = ['n', 'k', 'cycles', '1-cycles', '2-cycles', '3-cycles', '4-cycles', '5-cycles',
                      'max_cycle_len']
        log_lens = csv.DictWriter(f_log_lens, fieldnames=fieldnames)
        log_lens.writeheader()

        cycles = markov_process(n, k, p_aa, p_bb, p_ab, a_type, b_type, edges, f_log, log_lens)
        draw_number_of_cycles(cycles)


if __name__ == '__main__':
    main()
