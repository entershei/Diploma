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
    # Number of cycles with length is key. Cycle len is a number of fragile edges in it.
    cycles_m = {}
    # Number of cycles with different types
    cycles_with_edges_order = {}
    a_in_non_trivial_cycles = 0
    b_in_non_trivial_cycles = 0
    a_in_non_trivial_cycles_part = 0
    b_in_non_trivial_cycles_part = 0
    d = 0
    b = 0

    def __init__(
        self,
        cycles_m,
        cycles_with_edges_order,
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
        d,
        b,
    ):
        self.cycles_m = cycles_m
        self.cycles_with_edges_order = cycles_with_edges_order
        self.a_in_non_trivial_cycles = a_in_non_trivial_cycles
        self.b_in_non_trivial_cycles = b_in_non_trivial_cycles
        self.a_in_non_trivial_cycles_part = a_in_non_trivial_cycles_part
        self.b_in_non_trivial_cycles_part = b_in_non_trivial_cycles_part
        self.d = d
        self.b = b


def parse_logs_row(
    row,
    cycles_representatives,
    max_interesting_cycles_len,
    is_int,
):
    def to_type(value):
        return int(value) if is_int else float(value)

    cycles_m = {}
    cycles_with_edges_order = {}

    for cycle_type in cycles_representatives:
        cycles_with_edges_order[cycle_type] = to_type(row[cycle_type])

    for cycle_len in range(1, max_interesting_cycles_len):
        cycles_m[str(cycle_len)] = to_type(row[str(cycle_len)])

    a_in_non_trivial_cycles_part, b_in_non_trivial_cycles_part = -1, -1
    a_in_non_trivial_cycles = to_type(row["a_in_non_trivial_cycles"])
    b_in_non_trivial_cycles = to_type(row["b_in_non_trivial_cycles"])
    if "a_in_non_trivial_cycles_part" in row:
        a_in_non_trivial_cycles_part = to_type(row["a_in_non_trivial_cycles_part"])
        b_in_non_trivial_cycles_part = to_type(row["b_in_non_trivial_cycles_part"])
    d, b = -1, -1
    if "d" in row:
        d = to_type((row["d"]))
        b = to_type((row["b"]))

    return CyclesInfo(
        cycles_m,
        cycles_with_edges_order,
        a_in_non_trivial_cycles,
        b_in_non_trivial_cycles,
        a_in_non_trivial_cycles_part,
        b_in_non_trivial_cycles_part,
        d,
        b,
    )


# [min_len, max_len)
def generate_cycle_types_representative(min_len, max_len):
    to_representative = {}
    representatives = []
    for cycle_len in range(min_len, max_len):
        representatives += generate_cycle_types(min_len, max_len)
        all_cycles = generate_all_cycles(1, cycle_len, ["A", "B"])
        for cycle in all_cycles:
            to_representative[cycle] = "".join(["A"] * cycle.count("A") + ["B"] * cycle.count("B"))
    return to_representative, representatives


def generate_all_cycles(cur_len, max_cycle_len, cycles):
    if cur_len == max_cycle_len:
        return cycles
    cycles1 = list(map(lambda c: c + "A", cycles))
    cycles2 = list(map(lambda c: c + "B", cycles))
    return generate_all_cycles(cur_len + 1, max_cycle_len, cycles1 + cycles2)


def define_cycles_representative(max_m):
    def for_m(m):
        all_cycles = generate_all_cycles(1, m, ["A", "B"])
        used = set()
        res = {}
        representatives_m = []

        for cycle in all_cycles:
            if cycle in used:
                continue
            representatives_m.append(cycle)
            for i in range(len(cycle)):
                shifted_cycle = cycle[i:] + cycle[:i]
                used.add(shifted_cycle)
                res[shifted_cycle] = cycle
            reversed_cycle = cycle[::-1]
            for i in range(len(cycle)):
                shifted_cycle = reversed_cycle[i:] + reversed_cycle[:i]
                used.add(shifted_cycle)
                res[shifted_cycle] = cycle
        return res, representatives_m

    result = {}
    representatives = []
    for cycle_len in range(1, max_m):
        representatives_for_types, cur_representatives = for_m(cycle_len)
        result.update(representatives_for_types)
        representatives += cur_representatives
    return result, representatives


def read_experiments_cycles_info(
    f_in,
    max_cycle_len_with_types,
    max_interesting_cycles_len,
    is_int,
    number_of_experiments=None,
    is_cycles_ordered=False,
):
    if is_cycles_ordered:
        _, cycles_representatives = define_cycles_representative(max_cycle_len_with_types)
    else:
        cycles_representatives = generate_cycle_types(1, max_cycle_len_with_types)
    experiments = []

    with open(f_in, "r", newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",", quotechar="|")
        experiment = []
        for row in reader:
            k = int(row["k"])
            cycle_info = parse_logs_row(
                row, cycles_representatives, max_interesting_cycles_len, is_int
            )

            if k == 0 and len(experiment) > 0:
                experiments.append(experiment)
                experiment = [cycle_info]
                if (
                    number_of_experiments is not None
                    and len(experiments) == number_of_experiments
                ):
                    break
            else:
                experiment.append(cycle_info)

        if number_of_experiments is None or len(experiments) < number_of_experiments:
            experiments.append(experiment)

        csvfile.close()

    return experiments


def log_experiments(
    experiments,
    file,
    open_mode,
    max_interesting_cycles_len,
    cycles_representatives,
):
    # print("start log")
    file_exists = os.path.isfile(file)

    with open(file, open_mode, newline="") as f_log:
        cycle_lens = list(map(lambda l: str(l), range(1, max_interesting_cycles_len)))

        fieldnames = (
            ["k"]
            + [
                "a_in_non_trivial_cycles",
                "b_in_non_trivial_cycles",
                "a_in_non_trivial_cycles_part",
                "b_in_non_trivial_cycles_part",
                "d",
                "b",
            ]
            + cycle_lens
            + cycles_representatives
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
                "d": info.d,
                "b": info.b,
            }
            for cycle_type in cycles_representatives:
                if cycle_type in info.cycles_with_edges_order:
                    cur_result[cycle_type] = info.cycles_with_edges_order[cycle_type]
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
