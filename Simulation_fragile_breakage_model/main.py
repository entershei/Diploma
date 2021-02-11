import csv
import random

import matplotlib.pyplot as plt
import numpy as np


def read_number_of_fragile_regions_and_steps():
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


def update_cycles(edge1, edge2, joining_type, edges, colors, new_color):
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
            return new_color
        else:
            set_color(x, edges, colors, new_color)
            return new_color + 1
    else:
        set_color(x, edges, colors, colors[x])
        return new_color


def update_edges(edge1, edge2, join_type):
    x, y = edge1
    u, v = edge2

    if join_type == 'xu':
        new_edges = [[x, u], [y, v]]
    else:
        new_edges = [[x, v], [y, u]]

    random.shuffle(new_edges)
    return new_edges


class CyclesInfo:
    num_all_cycles = 0
    # cycles with length from 1 to 5
    num_n_cycles = [0] * 5
    max_len = 0

    def __init__(self, num_all_cycles, num_n_cycles, max_len):
        self.num_all_cycles = num_all_cycles
        self.num_n_cycles = num_n_cycles
        self.max_len = max_len


# Возвращает массив: количество циклов; [число циклов длины 1-5]; длина максимального цикла].
# Длина цикла измеряется в количестве хрупких ребер в нем.
def compute_cycles_lens(colors):
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

    return CyclesInfo(len(cycles), lens[:6], max_len)


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

    cnt_aa = 0
    cnt_bb = 0
    cnt_ab = 0

    steps_cycles_info.append(compute_cycles_lens(colors))

    for i in range(k):
        dcj_type = np.random.choice(['aa', 'bb', 'ab'], 1, p=probabilities)[0]
        joining_type = np.random.choice(['xu', 'xv'], 1)[0]

        if dcj_type == 'aa':
            cnt_aa += 1
            edge_ind1, edge_ind2 = np.random.choice(len(a_type), 2, replace=False)
            edge1 = a_type[edge_ind1]
            edge2 = a_type[edge_ind2]

            new_color = update_cycles(edge1, edge2, joining_type, edges, colors, new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            a_type[edge_ind1], a_type[edge_ind2] = new_edges
        elif dcj_type == 'bb':
            cnt_bb += 1
            edge_ind1, edge_ind2 = np.random.choice(len(b_type), 2, replace=False)
            edge1 = b_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color = update_cycles(edge1, edge2, joining_type, edges, colors, new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            b_type[edge_ind1], b_type[edge_ind2] = new_edges
        else:
            cnt_ab += 1
            edge_ind1 = np.random.choice(len(a_type), 1)[0]
            edge_ind2 = np.random.choice(len(b_type), 1)[0]

            edge1 = a_type[edge_ind1]
            edge2 = b_type[edge_ind2]

            new_color = update_cycles(edge1, edge2, joining_type, edges, colors, new_color)

            new_edges = update_edges(edge1, edge2, joining_type)
            a_type[edge_ind1], b_type[edge_ind2] = new_edges

        steps_cycles_info.append(compute_cycles_lens(colors))

        assert len_a == len(a_type) and len_b == len(b_type)

    return steps_cycles_info


def aggregate_cycles_info(experiments):
    num_experiments = len(experiments)
    num_steps_of_markov_process = len(experiments[0])
    num_interesting_cycles = 6

    aggregated_cycles_info = []
    for i in range(num_steps_of_markov_process):
        average_all_cycles_num = 0
        average_len_cycles = [0] * num_interesting_cycles
        average_max_cycles_len = 0
        for j in range(num_experiments):
            average_all_cycles_num += experiments[j][i].num_all_cycles
            for k in range(num_interesting_cycles):
                average_len_cycles[k] += experiments[j][i].num_n_cycles[k]
            average_max_cycles_len += experiments[j][i].max_len

        average_all_cycles_num /= num_experiments
        for k in range(num_interesting_cycles):
            average_len_cycles[k] /= num_experiments
        average_max_cycles_len /= num_experiments

        aggregated_cycles_info.append(CyclesInfo(average_all_cycles_num, average_len_cycles, average_max_cycles_len))

    return aggregated_cycles_info


def draw_number_of_cycles(cycles, title):
    k = len(cycles)
    x = list(range(k))
    plt.plot(x, cycles)
    plt.title(title)
    plt.xlabel('Number of swaps')
    plt.ylabel('Cycles')
    plt.grid()
    plt.show()


def log_aggregated_results(n, aggregated_cycles_info, f):
    with open(f, 'w', newline='') as f_log_lens:
        fieldnames = ['n', 'k', 'all-cycles', '1-cycles', '2-cycles', '3-cycles', '4-cycles', '5-cycles',
                      'max_cycle_len']
        log_lens = csv.DictWriter(f_log_lens, fieldnames=fieldnames)
        log_lens.writeheader()

        for step, info in enumerate(aggregated_cycles_info):
            log_lens.writerow(
                {'n': n, 'k': step, 'all-cycles': info.num_all_cycles, '1-cycles': info.num_n_cycles[1],
                 '2-cycles': info.num_n_cycles[2], '3-cycles': info.num_n_cycles[3],
                 '4-cycles': info.num_n_cycles[4], '5-cycles': info.num_n_cycles[5], 'max_cycle_len': info.max_len})


def main():
    # n, k = read_number_of_regions_and_steps()
    # p_aa, p_bb, p_ab = read_probabilities()
    n, p_aa, p_bb = 1000, 0.5, 0.5
    p_ab = 1 - p_aa - p_bb
    k = n * 2
    p_a_type_edges = 0.8
    p_split_a_b = [p_a_type_edges, 1 - p_a_type_edges]

    experiments = []
    for i in range(10):
        a_type, b_type, edges = split_fragile_edges(n, p_split_a_b)
        for j in range(10):
            experiments.append(markov_process(n, k, p_aa, p_bb, p_ab, a_type.copy(), b_type.copy(), edges.copy()))

    aggregated_cycles_info = aggregate_cycles_info(experiments)

    log_aggregated_results(n, aggregated_cycles_info, 'logs/log_lens_100_n1000_paa0_5_pbb0_5_betta0_8.csv')

    draw_number_of_cycles(list(map(lambda info: info.num_all_cycles, aggregated_cycles_info)),
                          'Number of cycles depends of number of swaps')
    draw_number_of_cycles(list(map(lambda info: info.num_n_cycles[1], aggregated_cycles_info)),
                          'Number of 1-cycles depends of number of swaps')
    draw_number_of_cycles(list(map(lambda info: info.num_n_cycles[2], aggregated_cycles_info)),
                          'Number of 2-cycles depends of number of swaps')
    draw_number_of_cycles(list(map(lambda info: info.num_n_cycles[3], aggregated_cycles_info)),
                          'Number of 3-cycles depends of number of swaps')
    draw_number_of_cycles(list(map(lambda info: info.num_n_cycles[4], aggregated_cycles_info)),
                          'Number of 4-cycles depends of number of swaps')
    draw_number_of_cycles(list(map(lambda info: info.num_n_cycles[5], aggregated_cycles_info)),
                          'Number of 5-cycles depends of number of swaps')


if __name__ == '__main__':
    main()
