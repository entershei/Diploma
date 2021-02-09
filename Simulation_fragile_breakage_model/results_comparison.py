import csv
import math


def read_logs(f):
    num_c1 = []
    with open(f, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        i = 0
        for row in reader:
            n = int(row['n'])
            k = int(row['k'])
            assert i == k

            c1 = float(row['1-cycles'])
            num_c1.append(c1)
            i += 1

    return num_c1, n


def compute_analytical_c1(gamma, p_aa, p_bb, betta):
    return betta * math.exp(-1 * gamma * (1 + p_aa - p_bb) / betta) +\
           (1 - betta) * math.exp(-1 * gamma * (1 - p_aa + p_bb) / (1 - betta))


# returns 100 * |analytical_c1 - real_c1| / real_c1
# betta = t / n, t is number of A-type edges
def compute_error(gamma, real_c1s, n, p_aa, p_bb, betta):
    k = int(gamma * n)

    real_c1 = real_c1s[k] / n

    return 100 * abs(compute_analytical_c1(gamma, p_aa, p_bb, betta) - real_c1) / real_c1


def write_error(error_depends_on_gamma, f):
    with open(f, 'w', newline='') as log:
        fieldnames = ['gamma', 'error']
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for gamma in error_depends_on_gamma.keys():
            log_lens.writerow(
                {'gamma': gamma, 'error': error_depends_on_gamma[gamma]})


# |analytical_c1 - real_c1| / real_c1
def compare_c1_error():
    real_c1s, n = read_logs('logs/log_lens_100_n1000_paa0_4_pbb0_35_betta0_7.csv')

    # gamma = k / n
    error_depends_on_gamma = {}
    gamma = 0.1
    step = 0.1
    while gamma < 2.1:
        error_depends_on_gamma[gamma] = compute_error(gamma, real_c1s, n, 0.4, 0.35, 0.7)
        gamma += step

    write_error(error_depends_on_gamma, 'logs/depends_on_gamma_n1000_paa0_4_pbb0_35_betta0_7.csv')


if __name__ == '__main__':
    compare_c1_error()
