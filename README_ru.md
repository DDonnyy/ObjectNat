# ObjectNat

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://img.shields.io/pypi/v/objectnat.svg)](https://pypi.org/project/objectnat/)
[![CI](https://github.com/DDonnyy/ObjectNat/actions/workflows/ci_pipeline.yml/badge.svg)](https://github.com/DDonnyy/ObjecNat/actions/workflows/ci_pipeline.yml)
[![codecov](https://codecov.io/gh/DDonnyy/ObjectNat/graph/badge.svg?token=K6JFSJ02GU)](https://codecov.io/gh/DDonnyy/ObjectNat)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/MIT)

<p align="center">
<img src="https://github.com/user-attachments/assets/ed0f226e-1728-4659-9e21-b4d499e703cd" alt="logo" width="400">
</p>

#### **ObjectNat** — это библиотека с открытым исходным кодом, разработанная командой **IDU** для пространственного анализа.

## Функции и как использовать

1. **[Изохроны и транспортная доступность](./examples/isochrone_generator.ipynb)** — Изохроны представляют собой области,
достижимые из исходной точки за заданное время по транспортной сети.
Эта функция позволяет анализировать транспортную доступность с использованием графов пешеходного, автомобильного,
общественного транспорта или их комбинации.

   Библиотека поддерживает несколько методов построения изохрон:
   - **Базовые изохроны**: отображают одну зону, достижимую за заданное время.
   - **Шаговые изохроны**: делят зону доступности на интервалы времени (например, 3, 5, 10 минут).

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_ways_15_min.png" alt="isochrone_ways_15_min" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_radius_15_min.png" alt="isochrone_radius_15_min" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_3points_radius_8_min.png" alt="isochrone_3points_radius_8_min" width="300">
   </p>
   <p align="center">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_ways_15_min.png" alt="stepped_isochrone_ways_15_min" width="300">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_radius_15_min.png" alt="stepped_isochrone_radius_15_min" width="300">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_separate_15_min.png" alt="stepped_isochrone_separate_15_min" width="300">
   </p>

2. **[Зоны покрытия](./examples/graph_coverage.ipynb)** — Функция генерации **зон покрытия** от набора исходных точек
с использованием транспортной сети. Вычисляет область, достижимую из каждой точки по **времени в пути** или **дистанции**,
затем строит полигоны с помощью **диаграмм Вороного** и обрезает их по заданной границе, если она указана.

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_time_10min.png" alt="coverage_zones_time_10min" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_distance_600m.png" alt="coverage_zones_distance_600m" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_radius_distance_800m.png" alt="coverage_zones_distance_radius_voronoi" width="350">
   </p>
   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_separate.png" alt="stepped_coverage_zones_separate" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_voronoi.png" alt="stepped_coverage_zones_voronoi" width="350">
   </p>

3. **[Анализ обеспеченности сервисами](./examples/calculate_provision.ipynb)** — Функция оценки обеспеченности жилых зданий
и их населения услугами (например, школы, поликлиники), которые имеют ограниченную **вместимость**
и заданный **порог доступности** (в минутах или метрах). Функция моделирует **баланс спроса и предложения**,
оценивая, насколько хорошо услуги удовлетворяют потребности близлежащих зданий в пределах допустимого времени.

   Библиотека также поддерживает:
   - **Перерасчёт** текущих результатов при изменении порога времени.
   - **Обрезку** результатов анализа по заданной зоне (например, границе района).

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_initial.png" alt="service_provision_initial" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_recalculated.png" alt="service_provision_recalculated" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_clipped.png" alt="service_provision_clipped" width="300">
   </p>

4. **[Анализ видимости](./examples/visibility_analysis.ipynb)** — Функция оценки видимости от заданной точки или множества
точек до близлежащих зданий в пределах заданного радиуса. Применяется для оценки визуальной доступности в городской среде.
Также реализован модуль для расчёта **зоны охвата** по видимости с использованием плотной сетки наблюдателей (рекомендуется ~1000 точек с шагом 10–20 метров).
Точки можно сгенерировать по транспортной сети и распределить по её рёбрам.

   Модуль включает:
   - **Быстрый приближённый метод** для больших объёмов данных.
   - **Точный метод** для локального детального анализа.

   <p align="center">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/visibility_comparison_methods.png" alt="visibility_comparison_methods" height="250">
     <img src="https://github.com/user-attachments/assets/b5b0d4b3-a02f-4ade-8772-475703cd6435" alt="visibility-catchment-area" height="250">
   </p>

5. **[Моделирование шума](./examples/noise_simulation.ipynb)** — Симуляция распространения шума от источников с учётом **препятствий**,
**растительности** и **факторов окружающей среды**.

   🔗 **[Подробное описание в Wiki](https://github.com/DDonnyy/ObjectNat/wiki/Симуляция-шумового-распространения)**

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/noise_simulation_1point.png" alt="noise_simulation_1point" width="400">
   </p>

6. **[Кластеризация точек](./examples/point_clusterization.ipynb)** — Функция построения **кластерных полигонов** по множеству точек на основе:
   - Минимального **расстояния** между точками.
   - Минимального **числа точек** в кластере.

   Также функция может рассчитывать **соотношение типов услуг** в каждом кластере для пространственного анализа состава услуг.

   <p align="center">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/building_clusters.png" alt="building_clusters" width="400">
   </p>

## Городские графы

Для достижения оптимальной производительности функций пространственного анализа ObjectNat рекомендуется использовать городские графы,
полученные с помощью библиотеки [IduEdu](https://github.com/DDonnyy/IduEdu).
**IduEdu** — это библиотека на Python с открытым исходным кодом, предназначенная для построения и обработки
сложных городских сетей на основе данных OpenStreetMap. 

**IduEdu** можно установить с помощью ``pip``:
```
pip install IduEdu
```

## Установка

**ObjectNat** можно установить с помощью ``pip``:

```
pip install ObjectNat
```

### Изменения конфигурации

```python
from objectnat import config

config.change_logger_lvl('INFO')  # Чтобы отключить отладочные сообщения
config.set_enable_tqdm(False)  # Чтобы отключить прогресс-бары tqdm
```

## Контакты

- [НЦКР](https://actcognitive.org/) — Национальный центр когнитивных разработок
- [ИДУ](https://idu.itmo.ru/) — Институт дизайна и урбанистики
- [Наталья Чичкова](https://t.me/nancy_nat) — менеджер проекта
- [Данила Олейников (Donny)](https://t.me/ddonny_dd) — ведущий инженер-разработчик

## Публикации
