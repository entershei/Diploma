import csv
import math


def read_logs(f, n_cycles):
    num_c_n = []
    with open(f, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        i = 0
        for row in reader:
            n = int(row['n'])
            k = int(row['k'])
            assert i == k

            c1 = float(row[n_cycles])
            num_c_n.append(c1)
            i += 1

    return num_c_n, n


def compute_analytical_c1(gamma, p_aa, p_bb, beta):
    return beta * math.exp(-1 * gamma * (1 + p_aa - p_bb) / beta) + \
           (1 - beta) * math.exp(-1 * gamma * (1 - p_aa + p_bb) / (1 - beta))


def compute_analytical_c2(gamma, p_aa, p_bb, beta):
    p1 = p_aa * gamma * math.exp(-2 * gamma * (1 + p_aa - p_bb) / beta)
    p2 = p_bb * gamma * math.exp(-2 * gamma * (1 + p_bb - p_aa) / (1 - beta))
    p3 = (1 - p_aa - p_bb) * gamma * math.exp(
        -1 * gamma * (1 - 2 * beta * p_aa + 2 * beta * p_bb + p_aa - p_bb) / (beta * (1 - beta)))

    return p1 + p2 + p3


# returns 100 * |analytical_c_n - real_c_n| / real_c_n
# beta = t / n, t is number of A-type edges
def compute_error(gamma, real_c_ns, n, p_aa, p_bb, beta, estimation_func):
    k = int(gamma * n)
    real_c_n = real_c_ns[k] / n
    analytical_c_n = estimation_func(gamma, p_aa, p_bb, beta)
    error = abs(analytical_c_n - real_c_n) / real_c_n

    return 100 * error


def write_error(error_depends_on_gamma, f):
    with open(f, 'w', newline='') as log:
        fieldnames = ['gamma', 'error']
        log_lens = csv.DictWriter(log, fieldnames=fieldnames)
        log_lens.writeheader()

        for gamma in error_depends_on_gamma.keys():
            log_lens.writerow(
                {'gamma': gamma, 'error': error_depends_on_gamma[gamma]})


# |analytical_c_n - real_c_n| / real_c_n
def compare_c_n_error(f_in, f_out, p_aa, p_bb, beta, n_cycles, compute_analytical_c_n):
    real_c_n_s, n = read_logs(f_in, n_cycles)

    # gamma = k / n
    error_depends_on_gamma = {}
    gamma = 0.1
    step = 0.1
    while gamma < 2.1:
        error_depends_on_gamma[gamma] = compute_error(gamma, real_c_n_s, n, p_aa, p_bb, beta, compute_analytical_c_n)
        gamma += step

    write_error(error_depends_on_gamma, f_out)


if __name__ == '__main__':
    # compare_c_n_error('logs/log_lens_5000_n1000_paa0_5_pbb0_45_beta0_5.csv',
    #                   'logs/error/1cycles/depends_on_gamma_n1000_paa0_5_pbb0_45_beta0_5.csv',
    #                   0.5, 0.45, 0.5, '1-cycles', compute_analytical_c1)
    file_end = '_n1000_paa0_4_pbb0_35_beta0_7.csv'
    compare_c_n_error(f_in='logs/log_lens_100' + file_end,
                      f_out='logs/error/2cycles/depends_on_gamma' + file_end,
                      p_aa=0.4, p_bb=0.35, beta=0.7, n_cycles='2-cycles',
                      compute_analytical_c_n=compute_analytical_c2)
