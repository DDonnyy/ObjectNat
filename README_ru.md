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

   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/isochrones.html)

2. **[Зоны покрытия](./examples/graph_coverage.ipynb)** — Функция генерации **зон покрытия** от набора исходных точек
с использованием транспортной сети. Вычисляет область, достижимую из каждой точки по **времени в пути** или **дистанции**,
затем строит полигоны с помощью **диаграмм Вороного** и обрезает их по заданной границе, если она указана.

   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/coverage.html)

3. **[Анализ обеспеченности сервисами](./examples/calculate_provision.ipynb)** — Функция оценки обеспеченности жилых зданий
и их населения услугами (например, школы, поликлиники), которые имеют ограниченную **вместимость**
и заданный **порог доступности** (в минутах или метрах). Функция моделирует **баланс спроса и предложения**,
оценивая, насколько хорошо услуги удовлетворяют потребности близлежащих зданий в пределах допустимого времени.

   Библиотека также поддерживает:
   - **Перерасчёт** текущих результатов при изменении порога времени.
   - **Обрезку** результатов анализа по заданной зоне (например, границе района).

   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/provision.html)

4. **[Анализ видимости](./examples/visibility_analysis.ipynb)** — Функция оценки видимости от заданной точки или множества
точек до близлежащих зданий в пределах заданного радиуса. Применяется для оценки визуальной доступности в городской среде.
Также реализован модуль для расчёта **зоны охвата** по видимости с использованием плотной сетки наблюдателей (рекомендуется ~1000 точек с шагом 10–20 метров).
Точки можно сгенерировать по транспортной сети и распределить по её рёбрам.

   Модуль включает:
   - **Быстрый приближённый метод** для больших объёмов данных.
   - **Точный метод** для локального детального анализа.
   
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/visibility.html)

5. **[Моделирование шума](./examples/noise_simulation.ipynb)** — Симуляция распространения шума от источников с учётом **препятствий**,
**растительности** и **факторов окружающей среды**.

   🔗 **[Подробное описание в Wiki](https://github.com/DDonnyy/ObjectNat/wiki/Симуляция-шумового-распространения)**
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/noise.html)

6. **[Кластеризация точек](./examples/point_clusterization.ipynb)** — Функция построения **кластерных полигонов** по множеству точек на основе:
   - Минимального **расстояния** между точками.
   - Минимального **числа точек** в кластере.

   Также функция может рассчитывать **соотношение типов услуг** в каждом кластере для пространственного анализа состава услуг.

   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/clustering.html)



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
