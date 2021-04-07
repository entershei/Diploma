NUMBER_OF_FRAGILE_EDGES = 1000
DIFFERENT_FRAGILE_EDGES_SPLITS = 200
MARKOV_PROCESS_EXPERIMENTS_ON_ONE_SPLIT = 5

EXPERIMENTS_FIELD_NAME = (
    "n"
    + str(NUMBER_OF_FRAGILE_EDGES)
    + "/"
    + str(DIFFERENT_FRAGILE_EDGES_SPLITS)
    + "_"
    + str(MARKOV_PROCESS_EXPERIMENTS_ON_ONE_SPLIT)
    + "_experiments/"
)

PROBABILITIES_WITH_ALPHA = [
    ("paa0,2_pbb0,6_alpha0,9.csv", 0.2, 0.6, 0.9),
    ("paa0,4_pbb0,35_alpha0,7.csv", 0.4, 0.35, 0.7),
    ("paa0,5_pbb0,5_alpha0,5.csv", 0.5, 0.5, 0.5),
    ("paa0,5_pbb0,5_alpha0,8.csv", 0.5, 0.5, 0.8),
    ("paa0,5_pbb0,45_alpha0,5.csv", 0.5, 0.45, 0.5),
    ("paa0,45_pbb0,45_alpha0,5.csv", 0.45, 0.45, 0.5),
]
