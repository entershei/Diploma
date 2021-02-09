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
    gammas, errors = read_logs('logs/depends_on_gamma_n1000_paa0_4_pbb0_35_betta0_7.csv')
    draw_error(gammas, errors, 'Percentage error of number of 1-cycles depends on gamma,' +
               '\nn = 1000, paa = 0.4, p_bb = 0.35, β = 0.7')


if __name__ == '__main__':
    main()
