import pandas as pd


def read_compartments():
    file = "3d_fragile_breakage_model/data/a_b_compartments.xlsx"

    df_compartments = pd.read_excel(file)
    chr_index = df_compartments.columns.get_loc("chr")
    start_index = df_compartments.columns.get_loc("start")
    end_index = df_compartments.columns.get_loc("end")
    h1_index = df_compartments.columns.get_loc("H1")

    a_b_compartments = []

    for i in range(len(df_compartments)):
        a_b_compartments.append(
            {
                "chr": df_compartments.iat[i, chr_index],
                "start": df_compartments.iat[i, start_index],
                "end": df_compartments.iat[i, end_index],
                "compartment": df_compartments.iat[i, h1_index],
            }
        )

    cnt_a = 0
    cnt_b = 0
    cnt_na = 0
    for i in a_b_compartments:
        if i["compartment"] == "A":
            cnt_a += 1
        elif i["compartment"] == "B":
            cnt_b += 1
        else:
            cnt_na += 1

    print("Sum of A regions:", cnt_a, "Mb")
    print("Sum of B regions:", cnt_b, "Mb")
    print("Sum of NA regions:", cnt_na, "Mb")

    return a_b_compartments


def build_breakpoint_graph():
    a_b_compartments = read_compartments()
    orthology_blocks =



if __name__ == "__main__":
    build_breakpoint_graph()
