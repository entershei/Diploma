import csv

import matplotlib.pyplot as plt
import numpy as np


def read_logs(f, cycle_types, to_summ, name_for_summ, to_rename):
    xs = []
    errors = []
    with open(f, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            xs.append(float(row['x']))
            error = {}
            for cycle_type in cycle_types:
                error[cycle_type] = float(row[cycle_type])
            for r in to_rename:
                error[r[1]] = error[r[0]]

            summ = 0
            for cycle_type in to_summ:
                summ += error[cycle_type]
            error[name_for_summ] = summ

            errors.append(error)

    return xs, errors


def draw_error(xs, errors, title, save_as):
    plt.plot(xs, errors)
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('Relative error')
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def draw_relative_errors(file_end, experiments, parameters_for_plot_name):
    # For 1-cycles
    cycle_types = ['A-cycles', 'B-cycles', '1-cycles']
    save_path = 'plots/relative_error/' + experiments + file_end + '/'
    xs, errors = read_logs('logs/relative_error/' + experiments + '1cycles/' +
                           'depends_on_x_' + file_end + '.csv', ['A-cycles', 'B-cycles', 'all'], [], '',
                           [['all', '1-cycles']])

    for cycle_type in cycle_types:
        draw_error(xs, list(map(lambda error: error[cycle_type], errors)),
                   'Relative error of number of ' + cycle_type + ' depends on x,\n' + parameters_for_plot_name,
                   save_path + cycle_type + '.png')

    # For 2-cycles
    cycle_types = ['AA-cycles', 'AB-cycles', 'BB-cycles', '2-cycles']
    save_path = 'plots/relative_error/' + experiments + file_end + '/'
    xs, errors = read_logs('logs/relative_error/' + experiments + '2cycles/' +
                           'depends_on_x_' + file_end + '.csv',
                           ['AA-cycles', 'AB-cycles', 'BB-cycles', 'all'], [], '',
                           [['all', '2-cycles']])

    for cycle_type in cycle_types:
        draw_error(xs, list(map(lambda error: error[cycle_type], errors)),
                   'Relative error of number of ' + cycle_type + ' depends on x,\n' + parameters_for_plot_name,
                   save_path + cycle_type + '.png')


def draw_number_of_cycles(cycles, title, save_as):
    k = len(cycles)
    x = list(range(k))
    plt.plot(x, cycles)
    plt.title(title)
    plt.xlabel('Number of swaps')
    plt.ylabel('Cycles')
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def draw_two_number_of_cycles(xs, real_cycles, analytical_cycles, title, save_as):
    plt.plot(xs, real_cycles, color='green', label='real')
    plt.plot(xs, analytical_cycles, color='blue', label='analytical')
    plt.legend(["Real", "Analytical"])
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('Cycles')
    plt.grid()
    # plt.show()
    plt.savefig(save_as)
    plt.close()


def read_cycles_info_logs(f):
    num_all_cycles = []
    different_cycles = {'A': [], 'B': [], 'AA': [], 'AB': [], 'BB': [], 'AAA': [], 'AAB': [], 'ABB': [], 'BBB': [],
                        'AAAA': [], 'AAAB': [], 'AABB': [], 'ABBB': [], 'BBBB': [],
                        'AAAAA': [], 'AAAAB': [], 'AAABB': [], 'AABBB': [], 'ABBBB': [], 'BBBBB': []}

    with open(f, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
        for row in reader:
            num_all_cycles.append(float(row['all-cycles']))
            different_cycles['A'].append(float(row['A-cycles']))
            different_cycles['B'].append(float(row['B-cycles']))
            different_cycles['AA'].append(float(row['AA-cycles']))
            different_cycles['AB'].append(float(row['AB-cycles']))
            different_cycles['BB'].append(float(row['BB-cycles']))
            different_cycles['AAA'].append(float(row['AAA-cycles']))
            different_cycles['AAB'].append(float(row['AAB-cycles']))
            different_cycles['ABB'].append(float(row['ABB-cycles']))
            different_cycles['BBB'].append(float(row['BBB-cycles']))
            different_cycles['AAAA'].append(float(row['AAAA-cycles']))
            different_cycles['AAAB'].append(float(row['AAAB-cycles']))
            different_cycles['AABB'].append(float(row['AABB-cycles']))
            different_cycles['ABBB'].append(float(row['ABBB-cycles']))
            different_cycles['BBBB'].append(float(row['BBBB-cycles']))
            different_cycles['AAAAA'].append(float(row['AAAAA-cycles']))
            different_cycles['AAAAB'].append(float(row['AAAAB-cycles']))
            different_cycles['AAABB'].append(float(row['AAABB-cycles']))
            different_cycles['AABBB'].append(float(row['AABBB-cycles']))
            different_cycles['ABBBB'].append(float(row['ABBBB-cycles']))
            different_cycles['BBBBB'].append(float(row['BBBBB-cycles']))

    return {'num_all_cycles': num_all_cycles, 'different_cycles': different_cycles,
            '1-cycles': np.add(different_cycles['A'], different_cycles['B']),
            '2-cycles': np.add(np.add(different_cycles['AA'], different_cycles['AB']), different_cycles['BB']),
            '3-cycles': np.add(np.add(different_cycles['AAA'], different_cycles['AAB']),
                               np.add(different_cycles['ABB'], different_cycles['BBB']))}


def draw_average_cycles(file_end, experiments):
    parameters = experiments + file_end
    f = 'logs/cycles_info/' + parameters + '.csv'
    save_path = 'plots/aggregated_cycles/' + parameters + '/'

    cycles_info = read_cycles_info_logs(f)

    draw_number_of_cycles(cycles_info['num_all_cycles'],
                          'Average number of cycles depends of number of swaps', save_path + 'all_cycles.png')
    draw_number_of_cycles(cycles_info['1-cycles'],
                          'Average number of 1-cycles depends of number of swaps', save_path + '1_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['A'],
                          'Average number of A-cycles depends of number of swaps', save_path + '1_A_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['B'],
                          'Average number of B-cycles depends of number of swaps', save_path + '1_B_cycles.png')

    draw_number_of_cycles(cycles_info['2-cycles'],
                          'Average number of 2-cycles depends of number of swaps', save_path + '2_cycles')
    draw_number_of_cycles(cycles_info['different_cycles']['AA'],
                          'Average number of AA-cycles depends of number of swaps', save_path + '2_AA_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['AB'],
                          'Average number of AB-cycles depends of number of swaps', save_path + '2_AB_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['BB'],
                          'Average number of BB-cycles depends of number of swaps', save_path + '2_BB_cycles.png')

    draw_number_of_cycles(cycles_info['3-cycles'],
                          'Average number of 3-cycles depends of number of swaps', save_path + '3_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['AAA'],
                          'Average number of AAA-cycles depends of number of swaps', save_path + '3_AAA_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['AAB'],
                          'Average number of AAB-cycles depends of number of swaps', save_path + '3_AAB_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['ABB'],
                          'Average number of ABB-cycles depends of number of swaps', save_path + '3_ABB_cycles.png')
    draw_number_of_cycles(cycles_info['different_cycles']['BBB'],
                          'Average number of BBB-cycles depends of number of swaps', save_path + '3_BBB_cycles.png')


def interesting_cycles_info(n, xs, cycles_info):
    interesting_info = []
    for x in xs:
        interesting_info.append(cycles_info[int(x * n)] / n)
    return interesting_info


def draw_average_with_analytical_cycles(file_end, experiments):
    n = 1000
    real_parameters = experiments + file_end
    save_path = 'plots/to_compare_number_of_cycles/' + real_parameters + '/'

    f_real = 'logs/cycles_info/' + real_parameters + '.csv'

    cycles_info = read_cycles_info_logs(f_real)

    # For 1-cycles
    cycle_types = ['A-cycles', 'B-cycles']
    xs, analytical_cycles = read_logs('logs/analytical_cycles/n1000/c1/' + file_end + '.csv',
                                      cycle_types, cycle_types, '1-cycles', [])

    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['1-cycles']),
                              list(map(lambda cycle: cycle['1-cycles'], analytical_cycles)),
                              'Normalized number of 1-cycles depends of x', save_path + '1_cycles.png')
    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['different_cycles']['A']),
                              list(map(lambda cycle: cycle['A-cycles'], analytical_cycles)),
                              'Normalized number of A-cycles depends of x', save_path + '1_A_cycles.png')
    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['different_cycles']['B']),
                              list(map(lambda cycle: cycle['B-cycles'], analytical_cycles)),
                              'Normalized number of B-cycles depends of x', save_path + '1_B_cycles.png')

    # For 2-cycles
    cycle_types = ['AA-cycles', 'AB-cycles', 'BB-cycles']
    xs, analytical_cycles = read_logs('logs/analytical_cycles/n1000/c2/' + file_end + '.csv',
                                      cycle_types, cycle_types, '2-cycles', [])

    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['2-cycles']),
                              list(map(lambda cycle: cycle['2-cycles'], analytical_cycles)),
                              'Normalized number of 2-cycles depends of x', save_path + '2_cycles.png')
    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['different_cycles']['AA']),
                              list(map(lambda cycle: cycle['AA-cycles'], analytical_cycles)),
                              'Normalized number of AA-cycles depends of x', save_path + '2_AA_cycles.png')
    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['different_cycles']['AB']),
                              list(map(lambda cycle: cycle['AB-cycles'], analytical_cycles)),
                              'Normalized number of AB-cycles depends of x', save_path + '2_AB_cycles.png')
    draw_two_number_of_cycles(xs, interesting_cycles_info(n, xs, cycles_info['different_cycles']['BB']),
                              list(map(lambda cycle: cycle['BB-cycles'], analytical_cycles)),
                              'Normalized number of BB-cycles depends of x', save_path + '2_BB_cycles.png')


def main():
    experiments = 'n1000/20_10_experiments/'
    file_end = 'paa0_5_pbb0_45_alpha0_5'
    parameters_for_plot_name = 'n = 1000, paa = 0.5, p_bb = 0.45, α = 0.5'

    draw_relative_errors(file_end, experiments, parameters_for_plot_name)
    draw_average_cycles(file_end, experiments)
    draw_average_with_analytical_cycles(file_end, experiments)


if __name__ == '__main__':
    main()
