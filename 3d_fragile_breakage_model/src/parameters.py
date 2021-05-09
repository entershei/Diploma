NUMBER_OF_FRAGILE_EDGES = 1000
NUMBER_OF_EXPERIMENTS = 3000
NUMBER_OF_STEPS = 3000
PART = NUMBER_OF_FRAGILE_EDGES / 20 + 1

EXPERIMENTS_IN_ONE_BUNCH = 100

MAX_POSSIBLE_CYCLES_LEN = 900

# p_aa, p_bb, a_type_edges_proportion
PROBABILITIES_WITH_ALPHA = [
    ("paa0,57142857142_pbb0,28571428571_alpha0,5", 0.57142857142, 0.28571428571, 0.5),
    ("paa0,4_pbb0,4_alpha0,5", 0.4, 0.4, 0.5),
    (
        "paa0,57142857142_pbb0,28571428571_alpha0,33333333333",
        0.57142857142,
        0.28571428571,
        0.33333333333,
    ),
    ("paa0,4_pbb0,4_alpha0,33333333333", 0.4, 0.4, 0.33333333333),
    ("paa0,48_pbb0,45_alpha0,53", 0.48, 0.45, 0.53),
    ("paa0,333_pbb0,333_alpha0,5", 0.333, 0.333, 0.5),
    ("paa0,4_pbb0,35_alpha0,7", 0.4, 0.35, 0.7),
    ("paa0,5_pbb0,5_alpha0,5", 0.5, 0.5, 0.5),
    ("paa0,5_pbb0,5_alpha0,8", 0.5, 0.5, 0.8),
    ("paa0,5_pbb0,45_alpha0,5", 0.5, 0.45, 0.5),
    ("paa0,45_pbb0,45_alpha0,5", 0.45, 0.45, 0.5),
    ("paa0,5_pbb0,4_alpha0,45", 0.5, 0.4, 0.45),
    ("paa0,45_pbb0,45_alpha0,48", 0.45, 0.45, 0.48),
]
