def generate_cycle_types(min_len, max_len):
    cycles = []
    for cur_len in range(min_len, max_len + 1):
        for cnt_a in range(cur_len, -1, -1):
            cycles.append("A" * cnt_a + "B" * (cur_len - cnt_a))
    return cycles


def generate_cycle_names(min_len, max_len):
    return map(
        lambda cycle_type: cycle_type + "-cycles",
        generate_cycle_types(min_len, max_len),
    )
