# Diploma

[dcj.py](fragile_breakage_model/src/dcj.py) записывает результаты Марковского процесса в logs/experiments.
Из-за большого количества экспериментов, данные о них занимают гигабайты, поэтому архив с экспериментами
[выложен](https://drive.google.com/drive/folders/1uALp2yIs_-3Dpz0_tv1styecyQEbf8A4?usp=sharing) на Гугл диск.

[aggregate_cycles_info.py](fragile_breakage_model/src/aggregate_cycles_info.py) обрабатывает результаты экспериментов из
logs/experiments и записывает агрегированные результаты в [logs/cycles_info/](fragile_breakage_model/logs/cycles_info/).

В [build_breakpoint_graph/](build_breakpoint_graph/) используется код из репозитория
[github.com/ctlab/true-dist-infer](https://github.com/ctlab/true-dist-infer) для построения графа точек разрыва
(breakpoint graph) из реальных данных в inferCARs формате.