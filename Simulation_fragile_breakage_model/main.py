import numpy as np
import matplotlib.pyplot as plt


def read_number_of_regions_and_steps():
    print('Enter number of regions (should be an even number)')
    n = int(input())
    assert n % 2 == 0

    print('Enter number of steps')
    k = int(input())

    return n, k


def read_probabilities():
    print('Enter probability of rearrangement between fragile regions:')
    p_ff = float(input())

    print('Enter probability of rearrangement between solid regions:')
    p_ss = float(input())

    # Probability of rearrangement between fragile and solid regions.
    p_fs = 1 - p_ss - p_ff

    assert p_ff >= 0 and p_ss >= 0 and p_fs >= 0

    return p_ff, p_ss, p_fs


def update_cycles(a, b, positions, cnt_cycles, colors, new_color):
    if positions[a] == b:
        a, b = b, a

    if colors[a] == colors[b]:
        if positions[b] == a:
            positions[b] = positions[a]
            positions[a] = a
            colors[a] = new_color
        else:
            positions[a], positions[b] = positions[b], positions[a]

            colors[a] = new_color
            cur = positions[a]
            while cur != a:
                colors[cur] = new_color
                cur = positions[cur]
        return cnt_cycles + 1, new_color + 1
    else:
        positions[a], positions[b] = positions[b], positions[a]
        colors[b] = colors[a]
        cur = positions[b]
        while cur != b:
            colors[cur] = colors[b]
            cur = positions[cur]

        return cnt_cycles - 1, new_color


def markov_process(n, k, p_ff, p_ss, p_fs, f_log):
    positions = list(range(-1, n))
    probabilities = [p_ff, p_ss, p_fs]

    cycles = [n]
    colors = list(range(n))
    cnt_cycles = n
    v_out = list(range(n))
    new_color = n

    for i in range(k):
        p = np.random.choice(3, 1, p=probabilities)

        if p == 0:
            a = np.random.randint(1, n // 2)
            a = a * 2 - 1
            b = np.random.randint(1, n // 2)
            b = b * 2 - 1
            regions = 'ff'
        elif p == 1:
            a = np.random.randint(1, n // 2)
            a = a * 2
            b = np.random.randint(1, n // 2)
            b = b * 2
            regions = 'ss'
        else:
            a = np.random.randint(1, n // 2)
            a = a * 2 - 1
            b = np.random.randint(1, n // 2)
            b = b * 2
            regions = 'fs'

        cnt_cycles, new_color = update_cycles(a, b, v_out, cnt_cycles, colors, new_color)
        cycles.append(cnt_cycles)

        # f_log.write('Regions = ' + regions + ' ' + str(a) + ' ' + str(b) + '\n')
        # f_log.write('Number of cycles = ' + str(cnt_cycles) + '\n')

        # f_log.write('Positions:\n')
        # for position in positions:
        #     f_log.write(str(position) + ' ')
        # f_log.write('\n\n')

    return cycles


def draw_number_of_cycles(cycles):
    k = len(cycles)
    x = list(range(k))
    plt.plot(x, cycles)
    plt.title('Number of cycles depends of number of swaps')
    plt.xlabel('Number of swaps')
    plt.ylabel('Cycles')
    plt.show()


def main():
    # n, k = read_number_of_regions_and_steps()
    # p_ff, p_ss, p_fs = read_probabilities()
    n, k, p_ff, p_ss, p_fs = 10000, 100000, 0.8, 0.15, 0.05

    f_log = open('log', 'w')
    cycles = markov_process(n, k, p_ff, p_ss, p_fs, f_log)
    draw_number_of_cycles(cycles)


if __name__ == '__main__':
    main()
