# ObjectNat 

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://img.shields.io/pypi/v/objectnat.svg)](https://pypi.org/project/objectnat/)
[![CI](https://github.com/DDonnyy/ObjecNat/actions/workflows/ci_pipeline.yml/badge.svg)](https://github.com/DDonnyy/ObjecNat/actions/workflows/ci_pipeline.yml)
[![Coverage](https://codecov.io/gh/DDonnyy/ObjecNat/graph/badge.svg?token=VN8CBP8ZW3)](https://codecov.io/gh/DDonnyy/ObjecNat)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/MIT)

<p align="center">
<img src="https://github.com/user-attachments/assets/d3878cce-8eba-4f96-8458-9a798d436120" alt="logo" width="400">
</p>

#### **ObjectNat** — это библиотека с открытым исходным кодом, предназначенная для геопространственного анализа, созданная командой **IDU**

## Компоненты ObjectNat

- [IduEdu](https://github.com/DDonnyy/IduEdu): `IduEdu` предоставляет функции для работы с графами
- [population-restorator](https://github.com/kanootoko/population-restorator): `restorator` предоставляет функции по расселению города

## Функции и как использовать

1. **[Граф города из OSM (IduEdu)](./examples/get_any_graph.ipynb)** — Функции для сборки графа дорог, пешеходных путей и общественного транспорта
из OpenStreetMap (OSM) и создания интермодального графа.

   <img src="https://github.com/user-attachments/assets/8dc98da9-8462-415e-8cc8-bdfca788e206" alt="IntermodalGraph" height="250">

2. **[Матрица смежности](./examples/calculate_adjacency_matrix.ipynb)** — Расчет матрицы смежности на основе переданного графа и типа весов рёбер
(время или расстояние). Графы можно получить, используя предыдущий пример.

3. **[Изохроны, транспортная доступность](./examples/isochrone_generator.ipynb)** — Функция для генерации изохрон для
анализа транспортной доступности от заданных начальных координат. Изохроны можно построить на основе графов
пешеходного, автомобильного или общественного транспорта.

   <img src="https://github.com/user-attachments/assets/37f308a5-db56-497d-b080-4edef3584fe5" alt="isochrones" height="250">

4. **[Восстановление населения](./examples/restore_population.ipynb)** — Функция для расселения населения в переданный слой жилых зданий.
Эта функция распределяет людей по домам на основе общей численности населения города и жилой площади каждого дома.

5. **[Обеспеченность услугами](./examples/calculate_provision.ipynb)** — Функция для расчета обеспеченности жилых зданий и населения
услугами на основе гравитационной модели.

   <img src="https://github.com/user-attachments/assets/5f2b3c55-9a02-4d70-80f4-503b77023eda" alt="ProvisionSchools" height="250">

6. **[Анализ видимости](./examples/visibility_analysis.ipynb)** — Функция для быстрого расчета видимости с заданной точки(точек) до зданий на заданном 
расстоянии. Также в примере указан калькулятор зоны охвата видимости для больших городских территорий. 
Эта функция предназначена для работы с как минимум 1000 точками, расположенными на расстоянии 10-20 метров друг от друга
для оптимальных результатов. Точки могут быть сгенерированы с использованием дорожного графа и случайного распределения
точек вдоль рёбер.

   <img src="https://github.com/user-attachments/assets/2927ac86-01e8-4b0e-9ea8-72ad81c13cf5" alt="visibility-from-point" height="250">
   
   <img src="https://github.com/user-attachments/assets/b5b0d4b3-a02f-4ade-8772-475703cd6435" alt="visibility-catchment-area" height="250">

7. **[Симуляция шума](./examples/noise_simulation.ipynb)** - Функция для симуляции шумового распространения от
источник-а(ов) с учётом характеристик источника, препятствий и зелёных насаждений 
**[Больше информации на странице Wiki](https://github.com/DDonnyy/ObjectNat/wiki/Симуляция-шумового-распространения)**

   <img src="https://github.com/user-attachments/assets/dd185867-67c4-4d03-8905-d06dd1d36fb3" alt="noise_sim" height="250">

8. **[Кластеризация точек](./examples/point_clusterization.ipynb)** — Функция для генерации полигонов кластеров для заданных точек на основе
минимального расстояния и минимального количества точек на кластер. Опционально можно рассчитать относительное
соотношение между типами услуг внутри кластеров.

    <img src="https://github.com/user-attachments/assets/2a9ad722-87d2-4954-9612-5ac3765aa824" alt="service-clusterization" height="250">

9. **[Жилые здания из OSM](./examples/download_buildings_from_osm.ipynb)** — Эта функция загружает геометрии зданий из OpenStreetMap (OSM) для указанной
территории и присваивает атрибуты каждому зданию. В частности, она определяет, является ли здание жилым
(атрибут `is_living`) и оценивает приблизительное количество жителей (атрибут`approximate_pop`).

   <img src="https://github.com/user-attachments/assets/d60dcd85-1a2e-4342-aae4-561aeda18858" alt="Living buildings" height="250">

## Установка

**ObjectNat** можно установить с помощью ``pip``:

```
pip install ObjectNat
```

### Изменения конфигурации

```python
from objectnat import config

config.set_timeout(10)  # Таймаут для запросов к Overpass
config.change_logger_lvl('INFO')  # Отключение всех сообщений отладки
config.set_enable_tqdm(False)  # Отключение всех индикаторов прогресса tqdm
config.set_overpass_url('http://your.overpass-api.de/interpreter/URL')
```

## Контакты

- [NCCR](https://actcognitive.org/) — Национальный центр когнитивных исследований
- [IDU](https://idu.itmo.ru/) — Институт дизайна и урбанистики
- [Наталья Чичкова](https://t.me/nancy_nat) — менеджер проекта
- [Данила Олейников (Donny)](https://t.me/ddonny_dd) — ведущий инженер-разработчик

## Публикации
