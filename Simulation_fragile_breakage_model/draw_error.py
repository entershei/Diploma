import csv

import matplotlib.pyplot as plt


def read_logs(f):
    gammas = []
    errors = []
    with open(f, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            gammas.append(float(row['gamma']))
            errors.append(float(row['error']))
    return gammas, errors


def draw_error(gammas, errors, title):
    plt.plot(gammas, errors)
    plt.title(title)
    plt.xlabel('γ')
    plt.ylabel('Percentage error')
    plt.grid()
    plt.show()


def main():
    gammas, errors = read_logs('logs/error/2cycles/depends_on_gamma_n1000_paa0_5_pbb_0_5_beta0_5.csv')
    draw_error(gammas, errors, 'Percentage error of number of 2-cycles depends on gamma,' +
               '\nn = 1000, paa = 0.5, p_bb = 0.5, β = 0.5')


if __name__ == '__main__':
    main()
