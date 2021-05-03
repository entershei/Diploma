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

    return a_b_compartments


def build_breakpoint_graph():
    a_b_compartments = read_compartments()


if __name__ == "__main__":
    build_breakpoint_graph()
