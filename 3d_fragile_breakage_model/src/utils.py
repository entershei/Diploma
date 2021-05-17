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
    a_in_non_trivial_cycles = 0
    b_in_non_trivial_cycles = 0
    a_in_non_trivial_cycles_part = 0
    b_in_non_trivial_cycles_part = 0

    def __init__(
        self,
        cycle_types,
        cycles_m,
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
    ):
        self.cycle_types = cycle_types
        self.cycles_m = cycles_m
        self.a_in_non_trivial_cycles = a_in_non_trivial_cycles
        self.b_in_non_trivial_cycles = b_in_non_trivial_cycles
        self.a_in_non_trivial_cycles_part = a_in_non_trivial_cycles_part
        self.b_in_non_trivial_cycles_part = b_in_non_trivial_cycles_part


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

    a_in_non_trivial_cycles_part, b_in_non_trivial_cycles_part = -1, -1
    if is_int:
        a_in_non_trivial_cycles = int(row["a_in_non_trivial_cycles"])
        b_in_non_trivial_cycles = int(row["b_in_non_trivial_cycles"])
        if "a_in_non_trivial_cycles_part" in row:
            a_in_non_trivial_cycles_part = int(row["a_in_non_trivial_cycles_part"])
            b_in_non_trivial_cycles_part = int(row["b_in_non_trivial_cycles_part"])
    else:
        a_in_non_trivial_cycles = float(row["a_in_non_trivial_cycles"])
        b_in_non_trivial_cycles = float(row["b_in_non_trivial_cycles"])
        if "a_in_non_trivial_cycles_part" in row:
            a_in_non_trivial_cycles_part = float(row["a_in_non_trivial_cycles_part"])
            b_in_non_trivial_cycles_part = float(row["b_in_non_trivial_cycles_part"])

    return CyclesInfo(
        cycle_types,
        cycles_m,
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
    )


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
    # print("start log")
    file_exists = os.path.isfile(file)

    with open(file, open_mode, newline="") as f_log:
        cycle_types = generate_cycle_types(1, max_cycle_len_with_types)
        cycle_lens = list(map(lambda l: str(l), range(1, max_interesting_cycles_len)))

        fieldnames = (
            ["k"]
            + [
                "a_in_non_trivial_cycles",
                "b_in_non_trivial_cycles",
                "a_in_non_trivial_cycles_part",
                "b_in_non_trivial_cycles_part",
            ]
            + cycle_types
            + cycle_lens
        )

        log_cycles = csv.DictWriter(f_log, fieldnames=fieldnames)
        if open_mode == "a" and not file_exists:
            log_cycles.writeheader()
        elif open_mode == "w":
            log_cycles.writeheader()

        for step, info in enumerate(experiments):
            cur_result = {
                "k": step,
                "a_in_non_trivial_cycles": info.a_in_non_trivial_cycles,
                "b_in_non_trivial_cycles": info.b_in_non_trivial_cycles,
                "a_in_non_trivial_cycles_part": info.a_in_non_trivial_cycles_part,
                "b_in_non_trivial_cycles_part": info.b_in_non_trivial_cycles_part,
            }
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
    # print("finish log")


def log_dictionaries(dictionaries, f):
    with open(f, "w", newline="") as f_log:
        log = csv.DictWriter(f_log, fieldnames=dictionaries[0].keys())
        log.writeheader()

        for dictionary in dictionaries:
            log.writerow(dictionary)


def read_logs(f):
    logs = []

    with open(f, "r", newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        for row in reader:
            logs.append(row)

        csvfile.close()
    return logs
