import csv

import os.path


def generate_cycle_type(cycle_len, cnt_a):
    return "A" * cnt_a + "B" * (cycle_len - cnt_a)


def generate_cycle_types_for_len(cycle_len):
    types = []
    for cnt_a in range(cycle_len, -1, -1):
        types.append(generate_cycle_type(cycle_len, cnt_a))
    return types


# [min_len; max_len)
def generate_cycle_types(min_len, max_len):
    types = []
    for cur_len in range(min_len, max_len):
        types += generate_cycle_types_for_len(cur_len)
    return types


class CyclesInfo:
    # Number of cycles with different types
    cycle_types = {}
    # Number of cycles with length is key. Cycle len is a number of fragile edges in it.
    cycles_m = {}

    def __init__(self, cycle_types, cycles_m):
        self.cycle_types = cycle_types
        self.cycles_m = cycles_m


def parse_logs_row(row, possible_cycle_types, max_interesting_cycles_len, is_int):
    cycle_types = {}
    cycles_m = {}
    for cycle_type in possible_cycle_types:
        if is_int:
            cycle_types[cycle_type] = int(row[cycle_type])
        else:
            cycle_types[cycle_type] = float(row[cycle_type])

    for cycle_len in range(1, max_interesting_cycles_len):
        if is_int:
            cycles_m[str(cycle_len)] = int(row[str(cycle_len)])
        else:
            cycles_m[str(cycle_len)] = float(row[str(cycle_len)])

    return CyclesInfo(cycle_types, cycles_m)


def read_experiments_cycles_info(
    f_in, max_cycle_len_with_types, max_interesting_cycles_len, is_int
):
    possible_cycle_types = generate_cycle_types(1, max_cycle_len_with_types)
    experiments = []

    with open(f_in, "r", newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        experiment = []
        for row in reader:
            k = int(row["k"])
            cycle_info = parse_logs_row(
                row, possible_cycle_types, max_interesting_cycles_len, is_int
            )

            if k == 0 and len(experiment) > 0:
                experiments.append(experiment)
                experiment = [cycle_info]
            else:
                experiment.append(cycle_info)

        experiments.append(experiment)

        csvfile.close()

    return experiments


def log_experiments(
    experiments, file, open_mode, max_cycle_len_with_types, max_interesting_cycles_len
):
    print("start log")
    file_exists = os.path.isfile(file)

    with open(file, open_mode, newline="") as f_log:
        cycle_types = generate_cycle_types(1, max_cycle_len_with_types)
        cycle_lens = list(map(lambda l: str(l), range(1, max_interesting_cycles_len)))

        fieldnames = ["k"] + cycle_types + cycle_lens

        log_cycles = csv.DictWriter(f_log, fieldnames=fieldnames)
        if open_mode == "a" and not file_exists:
            log_cycles.writeheader()
        elif open_mode == "w":
            log_cycles.writeheader()

        for step, info in enumerate(experiments):
            cur_result = {"k": step}
            for cycle_type in cycle_types:
                if cycle_type in info.cycle_types:
                    cur_result[cycle_type] = info.cycle_types[cycle_type]
                else:
                    cur_result[cycle_type] = 0

            for c_len in cycle_lens:
                if c_len in info.cycles_m:
                    cur_result[c_len] = info.cycles_m[c_len]
                else:
                    cur_result[c_len] = 0
            log_cycles.writerow(cur_result)
        f_log.close()
    print("finish log")
