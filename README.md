# Diploma

[simulation.py](3d_fragile_breakage_model/src/simulation.py) записывает результаты Марковского процесса в
[logs/experiments_100](3d_fragile_breakage_model/logs/experiments_100/), суммируя данные по 100 экспериментам.

[aggregate_cycles_info.py](3d_fragile_breakage_model/src/aggregate_cycles_info.py) обрабатывает результаты экспериментов
из [logs/experiments_100](3d_fragile_breakage_model/logs/experiments_100/) и записывает агрегированные результаты в
[logs/cycles_info/](3d_fragile_breakage_model/logs/cycles_info/).

В [build_breakpoint_graph/](build_breakpoint_graph/) используется код из репозитория
[github.com/ctlab/true-dist-infer](https://github.com/ctlab/true-dist-infer) для построения графа точек разрыва
(breakpoint graph) из реальных данных в inferCARs формате.

Данные о разбиении генома человека на [A/B compartments](3d_fragile_breakage_model/data/a_b_compartments.xlsx) взяты из
статьи [A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome / A. D. Schmitt [et al.] // Cell Reports. — 2016. — Volume 17, Issue 8, p. 2042-2059](https://doi.org/10.1016/j.celrep.2016.10.061).