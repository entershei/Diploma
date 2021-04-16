import csv

from results_comparison import compute_analytical_cycles_m
from utils import (
    generate_cycle_types)


def read_cycles_info(f, min_len, max_len):
    cycles_types = generate_cycle_types(min_len, max_len)
    num_c_n = []
    with open(f, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        i = 0
        for row in reader:
            n = int(row["n"])
            k = int(row["k"])
            assert i == k

            c_n = {}
            for cycle_type in cycles_types:
                c_n[cycle_type] = float(row[cycle_type])
            num_c_n.append(c_n)
            i += 1

    return num_c_n, n


def check_d_b(f):
    return 0


def compute_min_dist(x, p_aa, p_bb, alpha, eps):
    dist = 0.0
    m = 2
    delta = 100
    while delta > eps:
        old_res = dist
        dist += (m - 1) * compute_analytical_cycles_m(m, x, p_aa, p_bb, alpha)
        delta = dist - old_res
        # print(m, dist, delta)
        m += 1

    return dist, m


def main():
    p = list(map(lambda a: a / 10, range(1, 10)))
    for p_aa in p:
        for p_bb in p:
            if p_aa + p_bb < 1:
                for x in p:
                    for alpha in p:
                        d, iterations = compute_min_dist(x, p_aa, p_bb, alpha, eps=1e-4 * 4)
                        print(d, iterations, p_aa, p_bb, x, alpha)

    # check_d_b(get_parameters_as_string() + "paa0,5_pbb0,45_alpha0,5_10000_experiments.csv")


if __name__ == "__main__":
    main()
