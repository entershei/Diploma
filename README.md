# Diploma

[simulation.py](3d_fragile_breakage_model/src/simulation.py) записывает результаты Марковского процесса в
[logs/experiments_100](3d_fragile_breakage_model/logs/experiments_100/), суммируя данные по 100 экспериментам.

Данные, записанные о каждом эксперименте отдельно,
[выложены](https://drive.google.com/drive/folders/1uALp2yIs_-3Dpz0_tv1styecyQEbf8A4?usp=sharing) на Гугл диск из-за
большого объема данных. Так же на Гугл диск выложены данные некоторых экспериментов из
[logs/experiments_100](3d_fragile_breakage_model/logs/experiments_100/), которые занимают больше 100MB.

[aggregate_cycles_info.py](3d_fragile_breakage_model/src/aggregate_cycles_info.py) обрабатывает результаты экспериментов
из logs/experiments и записывает агрегированные результаты в
[logs/cycles_info/](3d_fragile_breakage_model/logs/cycles_info/).

В [build_breakpoint_graph/](build_breakpoint_graph/) используется код из репозитория
[github.com/ctlab/true-dist-infer](https://github.com/ctlab/true-dist-infer) для построения графа точек разрыва
(breakpoint graph) из реальных данных в inferCARs формате.