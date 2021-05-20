import random
from collections import defaultdict
import networkx as nx
import pandas as pd

from build_breakpoint_graph.src.utils.parsers import parse_to_df
from compute_statistics import compute_d_by_cycles, compute_b_by_cycles
from src.graphs.real_data_graph import (
    RealDataGraph,
    nodes_edges_one_sp,
    sp_blocks_from_df,
    uniq_predicate,
)
from true_evolutionary_distance import (
    find_true_evolution_dist_and_find_parameters1,
)
from utils import CyclesInfo, log_dictionaries, define_cycles_representative


def read_compartments(file):
    df_compartments = pd.read_excel(file)
    chr_index = df_compartments.columns.get_loc("chr")
    start_index = df_compartments.columns.get_loc("start")
    end_index = df_compartments.columns.get_loc("end")
    h1_index = df_compartments.columns.get_loc("H1")

    a_b_compartments = {}

    cnt_a = 0
    cnt_b = 0
    cnt_na = 0

    for i in range(len(df_compartments)):
        if df_compartments.iat[i, chr_index] not in a_b_compartments:
            a_b_compartments[df_compartments.iat[i, chr_index]] = []
        a_b_compartments[df_compartments.iat[i, chr_index]].append(
            {
                "start": df_compartments.iat[i, start_index],
                "end": df_compartments.iat[i, end_index],
                "compartment": df_compartments.iat[i, h1_index],
            }
        )
        if df_compartments.iat[i, h1_index] == "A":
            cnt_a += 1
        elif df_compartments.iat[i, h1_index] == "B":
            cnt_b += 1
        else:
            cnt_na += 1

    # print("Sum of A regions in compartments file:", cnt_a, "Mb")
    # print("Sum of B regions in compartments file:", cnt_b, "Mb")
    # print("Sum of NA regions in compartments file:", cnt_na, "Mb")
    # print()

    return a_b_compartments


def read_files_path():
    print("Enter a path to a file that contains A and B compartments:")
    compartments_file = input()

    print("Enter a path to a file with orthology blocks:")
    blocks_file = input()

    return compartments_file, blocks_file


def read_species(orthology_blocks):
    species = orthology_blocks.species.unique()
    # print("Choose two species from the list:")
    # must match the one for which the blocking was introduced
    print("Choose two species from the list to compare with homo_sapiens)")
    print(*species, sep=", ")

    print(
        "Format: species1, species2\n",
        "The second species must coincide with the one for which the partitioning into A/B compartments was",
        " introduced:",
    )

    species1, species2 = input().split()

    if species1 not in species:
        print("First species is not from the list")
        exit()

    if species2 not in species:
        print("Second species is not from the list")
        exit()
    return species1, species2


def add_types_to_fragile_edges(edges, a_b_compartments):
    def chr_to_num(chromosome):
        return int(chromosome) if chromosome != "X" else 23

    def find_chrs_ends():
        ends = {}
        for chromosome in a_b_compartments.keys():
            ends[chromosome] = a_b_compartments[chromosome][-1]["end"]
        return ends

    def compartments_for_edge(edge, cur_chr):
        cur_it_comp = it_compartment[cur_chr]

        cur_a_compartments = 0
        cur_b_compartments = 0
        cur_edge_inside_compartment = False
        while a_b_compartments[cur_chr][cur_it_comp]["end"] < edge["chr_beg"]:
            cur_it_comp += 1

        if (
            a_b_compartments[cur_chr][cur_it_comp]["start"] <= edge["chr_beg"]
            and edge["chr_end"] <= a_b_compartments[cur_chr][cur_it_comp]["end"] + 1
        ):
            cur_edge_inside_compartment = True

        while (
            cur_it_comp < len(a_b_compartments[cur_chr])
            and a_b_compartments[cur_chr][cur_it_comp]["start"] < edge["chr_end"]
        ):
            to_add = max(
                0,
                min(edge["chr_end"], a_b_compartments[cur_chr][cur_it_comp]["end"] + 1)
                - max(edge["chr_beg"], a_b_compartments[cur_chr][cur_it_comp]["start"]),
            )
            if a_b_compartments[cur_chr][cur_it_comp]["compartment"] == "A":
                cur_a_compartments += to_add
            elif a_b_compartments[cur_chr][cur_it_comp]["compartment"] == "B":
                cur_b_compartments += to_add

            cur_it_comp += 1
        cur_it_comp -= 1

        it_compartment[cur_chr] = cur_it_comp
        return cur_a_compartments, cur_b_compartments, cur_edge_inside_compartment

    it_compartment = defaultdict(lambda: 0)
    chrs_end = find_chrs_ends()

    cnt_edge_inside_compartment = 0
    cnt_edge_contains_multi_compartments = 0
    a_and_b = 0

    for cur_edge in edges:
        cur_chr = chr_to_num(cur_edge["chr"])

        if cur_edge["chr_beg"] > cur_edge["chr_end"]:
            a_compartments1, b_compartments1, _ = compartments_for_edge(
                {"chr_beg": cur_edge["chr_beg"], "chr_end": chrs_end[cur_chr]}, cur_chr
            )
            it_compartment[cur_chr] = 0
            a_compartments2, b_compartments2, _ = compartments_for_edge(
                {"chr_beg": 1, "chr_end": cur_edge["chr_end"]}, cur_chr
            )
            a_compartments = a_compartments1 + a_compartments2
            b_compartments = b_compartments1 + b_compartments2
        else:
            a_compartments, b_compartments, inside = compartments_for_edge(
                cur_edge, cur_chr
            )
            if inside:
                cnt_edge_inside_compartment += 1

        cur_edge["label"] = "A" if a_compartments > b_compartments else "B"
        if a_compartments == b_compartments:
            cur_edge["label"] = random.choice(["A", "B"])

        if a_compartments > 0 and b_compartments > 0:
            a_and_b += 1

        if a_compartments + b_compartments > 1e6:
            cnt_edge_contains_multi_compartments += 1

    # print("edge inside", cnt_edge_inside_compartment)
    # print("edge contains", cnt_edge_contains_multi_compartments)
    # print("A and B", a_and_b)
    #
    # cnt_a = defaultdict(lambda: 0)
    # cnt_b = defaultdict(lambda: 0)
    # for edge in edges:
    #     if edge["label"] == "A":
    #         cnt_a[edge["chr"]] += 1
    #     else:
    #         cnt_b[edge["chr"]] += 1

    # for chr_ in cnt_a.keys():
    #     print(chr_, cnt_a[chr_], cnt_b[chr_])


def len_of_blocks(orthology_blocks, species):
    len_blocks = 0
    species_blocks = orthology_blocks.loc[orthology_blocks["species"] == species]
    for _, block in species_blocks.iterrows():
        len_blocks += block.chr_end - block.chr_beg
    # print("Length of orthology blocks for", species, ":", len_blocks, "\n")


def get_compartments_and_orthology_blocks(compartments_file, blocks_file):
    a_b_compartments = read_compartments(compartments_file)
    orthology_blocks = parse_to_df(blocks_file)
    return a_b_compartments, orthology_blocks


def build_breakpoint_graph(species1, species2, a_b_compartments, orthology_blocks):
    print("Species:", species1 + ",", species2)
    print()

    len_of_blocks(orthology_blocks, species1)

    allowed_blocks = sp_blocks_from_df(orthology_blocks, species1)
    allowed_blocks = list(
        filter(
            uniq_predicate((sp_blocks_from_df(orthology_blocks, species2))),
            allowed_blocks,
        )
    )

    nodes1, edges1 = nodes_edges_one_sp(
        "black", species1, orthology_blocks, allowed_blocks, True
    )
    nodes2, edges2 = nodes_edges_one_sp(
        "red", species2, orthology_blocks, allowed_blocks, True
    )
    add_types_to_fragile_edges(edges2, a_b_compartments)

    g = RealDataGraph()
    g.add_nodes_and_edges(nodes1 + nodes2, edges1 + edges2)

    return g


def get_graph_statistic(g):
    def to_order(edges):
        def dfs(v, edge_to):
            used[v] = True
            ordered_edges.append(edge_to)
            for u_edge in adjacency[v]:
                if not used[u_edge[1]]:
                    dfs(u_edge[1], u_edge)

        adjacency = {}
        used = {}
        for i, edge in enumerate(edges):
            used[edge[0]] = False
            used[edge[1]] = False
            if edge[0] not in adjacency:
                adjacency[edge[0]] = [edge]
            else:
                adjacency[edge[0]].append(edge)

            if edge[1] not in adjacency:
                adjacency[edge[1]] = [(edge[1], edge[0], edge[2])]
            else:
                adjacency[edge[1]].append((edge[1], edge[0], edge[2]))

        used[edges[0][0]] = True
        ordered_edges = []
        dfs(edges[0][1], edges[0])

        for connected_to_first in adjacency[edges[0][0]]:
            if connected_to_first[1] != edges[0][1] or connected_to_first[2] != edges[0][2]:
                ordered_edges.append((connected_to_first[1], connected_to_first[0], connected_to_first[2]))
                break

        return ordered_edges

    def to_edges_types(edges):
        edges = to_order(edges)
        edges_types = []
        for edge in edges:
            if edge[2] == "A" or edge[2] == "B":
                edges_types.append(edge[2])
        return edges_types

    def cycle_to_represent(edges):
        if len(edges) < max_cycles_len_with_order:
            return to_represent["".join(edges)]
        cycle_type = ["A"] * edges.count("A") + ["B"] * edges.count("B")
        return "".join(cycle_type)

    max_cycles_len_with_order = 6
    to_represent, _ = define_cycles_representative(max_cycles_len_with_order)
    statistic = {
        "cycles": defaultdict(lambda: 0),
        "cycles_with_types": defaultdict(lambda: 0),
        "total A-edges": 0,
        "total B-edges": 0,
        "A-edges in non-trivial cycles": 0,
        "B-edges in non-trivial cycles": 0,
    }
    for component in nx.connected_components(g):
        cycle = g.subgraph(component)
        if len(cycle.nodes) == len(cycle.edges):
            cycle_len = len(cycle.edges) // 2
            statistic["cycles"][str(cycle_len)] += 1
            all_edges = to_edges_types(list(cycle.edges(data="label")))
            statistic["cycles_with_types"][cycle_to_represent(all_edges)] += 1

            if cycle_len > 1:
                statistic["A-edges in non-trivial cycles"] += all_edges.count("A")
                statistic["B-edges in non-trivial cycles"] += all_edges.count("B")
            statistic["total A-edges"] += all_edges.count("A")
            statistic["total B-edges"] += all_edges.count("B")

    return statistic


def print_graph_statistic(graph_statistic):
    total_cycles = 0
    for cycle_len in graph_statistic["cycles"].keys():
        total_cycles += graph_statistic["cycles"][cycle_len]
    print("Total cycles:", total_cycles)

    for cycle_len in sorted(graph_statistic["cycles"].keys()):
        print("  " + str(cycle_len) + "-cycles:", graph_statistic["cycles"][cycle_len])

    print("Cycles with types:")
    for cycle_type in sorted(
        graph_statistic["cycles_with_types"].keys(), key=lambda c_type: len(c_type)
    ):
        c_type = cycle_type
        if len(c_type) > 10:
            c_type = (
                str(cycle_type.count("A")) + "-A," + str(cycle_type.count("B")) + "-B"
            )
        print(
            "  " + c_type + "-cycles:",
            graph_statistic["cycles_with_types"][cycle_type],
        )

    print()
    print(
        "Total A-edges:",
        graph_statistic["total A-edges"],
    )
    print(
        "Total B-edges:",
        graph_statistic["total B-edges"],
    )

    print(
        "A-edges in non-trivial cycles:",
        graph_statistic["A-edges in non-trivial cycles"],
    )
    print(
        "B-edges in non-trivial cycles:",
        graph_statistic["B-edges in non-trivial cycles"],
    )
    print()


def find_true_evolution_dist(graph_statistic):
    cycles_m = dict(graph_statistic["cycles"])
    graph_cycles_info = CyclesInfo(
        cycles_m,
        dict(graph_statistic["cycles_with_types"]),
        graph_statistic["A-edges in non-trivial cycles"],
        graph_statistic["B-edges in non-trivial cycles"],
        -1,
        -1,
        compute_d_by_cycles(cycles_m),
        compute_b_by_cycles(cycles_m),
    )
    sum_to_in_error = 10
    for i in range(2, sum_to_in_error):
        if str(i) not in graph_cycles_info.cycles_m:
            graph_cycles_info.cycles_m[str(i)] = 0

    max_cycles_len_with_order = 6
    to_represent, _ = define_cycles_representative(max_cycles_len_with_order)
    dist_info = find_true_evolution_dist_and_find_parameters1(
        graph_cycles_info, to_represent
    )

    print("Estimated true evolutionary distance:", dist_info["estimated_true_dist"])
    print(
        "Parameters at which the distance is reached:\n    p_aa:",
        dist_info["best_p_aa"],
        "p_ab:",
        1 - dist_info["best_p_aa"] - dist_info["best_p_bb"],
        "p_bb:",
        dist_info["best_p_bb"],
        "alpha:",
        dist_info["best_alpha"],
        "beta:",
        1 - dist_info["best_alpha"],
    )
    print("Min error:", dist_info["min_error"])
    print("Empirical min distance:", dist_info["empirical_min_dist"])
    print()
    return dist_info


def estimate_distance_for_many_species(
    species_pairs, a_b_compartments, orthology_blocks
):
    for species_pair in species_pairs:
        species1, species2 = species_pair
        g = build_breakpoint_graph(
            species1, species2, a_b_compartments, orthology_blocks
        )
        graph_statistic = get_graph_statistic(g)
        print_graph_statistic(graph_statistic)
        dist_info = find_true_evolution_dist(graph_statistic)
        dist_info["species1"] = species1
        dist_info["species2"] = species2

        log_dictionaries(
            [{**dist_info, **graph_statistic}],
            "3d_fragile_breakage_model/logs/true_evolution_distance_found_parameters/real_genomes/"
            + species1
            + "_"
            + species2
            + ".csv",
        )


def main():
    # compartments_file, blocks_file = read_files_path()
    compartments_file = "3d_fragile_breakage_model/data/a_b_compartments.xlsx"
    blocks_file = "3d_fragile_breakage_model/data/orthology_blocks.txt"

    a_b_compartments, orthology_blocks = get_compartments_and_orthology_blocks(
        compartments_file, blocks_file
    )
    # species1, species2 = read_species(orthology_blocks)
    estimate_distance_for_many_species(
        [
            # ("macaca_mulatta", "homo_sapiens"),
            ("mus_musculus", "homo_sapiens"),
            # ("rattus_norvegicus", "homo_sapiens"),
            # ("monodelphis_domestica", "homo_sapiens"),
        ],
        a_b_compartments,
        orthology_blocks,
    )


if __name__ == "__main__":
    main()
