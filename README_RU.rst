ObjectNat
=========

Object-oriented Network Analysis Tools
--------------------------------------

.. |badge-black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Стиль кода: black

.. |badge-pypi| image:: https://img.shields.io/pypi/v/objectnat.svg
   :target: https://pypi.org/project/objectnat/
   :alt: Версия PyPI

.. |badge-ci| image:: https://github.com/DDonnyy/ObjectNat/actions/workflows/ci_pipeline.yml/badge.svg
   :target: https://github.com/DDonnyy/ObjectNat/actions/workflows/ci_pipeline.yml
   :alt: CI статус

.. |badge-codecov| image:: https://codecov.io/gh/DDonnyy/ObjectNat/graph/badge.svg?token=K6JFSJ02GU
   :target: https://codecov.io/gh/DDonnyy/ObjectNat
   :alt: Покрытие тестами

.. |badge-license| image:: https://img.shields.io/badge/license-BSD--3--Clause-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: Лицензия

.. |badge-docs| image:: https://img.shields.io/badge/docs-latest-4aa0d5?logo=readthedocs
   :target: https://iduclub.github.io/ObjectNat/
   :alt: Документация

|badge-black| |badge-pypi| |badge-ci| |badge-codecov| |badge-license| |badge-docs|

`README (English) <README.rst>`_

.. image:: https://raw.githubusercontent.com/IDUclub/ObjectNat/main/docs/_static/ONlogo.svg
   :align: center
   :width: 400
   :alt: ObjectNat logo

----

**ObjectNat** — это библиотека с открытым исходным кодом, разработанная командой **IDU**
для пространственного и сетевого анализа в городских исследованиях.
Библиотека предоставляет инструменты для анализа **доступности**, **видимости**,
**распространения шума** и **обеспеченности сервисами**.

----

Основные функции
----------------

Каждая функция сопровождается **примером в Jupyter Notebook** и **документацией**.

1. **Изохроны и транспортная доступность**

   Изохроны представляют собой области, достижимые из исходной точки за заданное время по транспортной сети.
   Эта функция позволяет анализировать транспортную доступность с использованием графов пешеходного, автомобильного,
   общественного транспорта или их комбинации.

   Библиотека поддерживает несколько методов построения изохрон:

   - **Базовые изохроны**: отображают одну зону, достижимую за заданное время.
   - **Шаговые изохроны**: делят зону доступности на интервалы времени (например, 3, 5, 10 минут).

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/isochrones.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/isochrones.html>`__

2. **Зоны покрытия**

   Функция генерации **зон покрытия** от набора исходных точек с использованием транспортной сети. Вычисляет область,
   достижимую из каждой точки по **времени в пути** или **дистанции**, затем строит полигоны с помощью
   **диаграмм Вороного** и обрезает их по заданной границе, если она указана.

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/coverage.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/coverage.html>`__

3. **Анализ обеспеченности сервисами**

   Функция оценки обеспеченности жилых зданий и их населения услугами (например, школы, поликлиники),
   которые имеют ограниченную **вместимость** и заданный **порог доступности** (в минутах или метрах).
   Функция моделирует **баланс спроса и предложения**, оценивая, насколько хорошо услуги удовлетворяют потребности
   близлежащих зданий в пределах допустимого времени.

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/provision.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/provision.html>`__

4. **Анализ видимости**

   Функция оценки видимости от заданной точки или множества точек до близлежащих зданий в пределах заданного радиуса.
   Применяется для оценки визуальной доступности в городской среде. Также реализован модуль для расчёта **зоны охвата**
   по видимости с использованием плотной сетки наблюдателей (рекомендуется ~1000 точек с шагом 10–20 метров).
   Точки можно сгенерировать по транспортной сети и распределить по её рёбрам.

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/visibility.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/visibility.html>`__

5. **Моделирование шума**

   Симуляция распространения шума от источников с учётом **препятствий**, **растительности** и **факторов окружающей среды**.

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/noise.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/noise.html>`__
   🧠 `Подробное описание <https://github.com/DDonnyy/ObjectNat/wiki/Noise-simulation>`__

6. **Кластеризация точек**

   Функция построения **кластерных полигонов** по множеству точек на основе:

   - Минимального **расстояния** между точками.
   - Минимального **числа точек** в кластере.

   Также функция может рассчитывать **соотношение типов услуг** в каждом кластере для пространственного анализа состава услуг.

   📘 `Пример <https://iduclub.github.io/ObjectNat/methods/examples/clustering.html>`__
   🔗 `Документация <https://iduclub.github.io/ObjectNat/methods/clustering.html>`__

----

Городские графы с помощью *IduEdu*
----------------------------------

Для оптимальной работы **ObjectNat** рекомендуется использовать графы,
созданные библиотекой `IduEdu <https://github.com/IDUclub/IduEdu>`__.

**IduEdu** — это библиотека на Python с открытым исходным кодом, предназначенная для построения и обработки
сложных городских сетей на основе данных OpenStreetMap.


**IduEdu** можно установить с помощью ``pip``::

    pip install IduEdu

Пример использования::

    from iduedu import get_4326_boundary, get_intermodal_graph

    poly = get_4326_boundary(osm_id=1114252)
    G_intermodal = get_intermodal_graph(territory=poly, clip_by_territory=True)

----

Установка
---------

**ObjectNat** можно установить с помощью ``pip``::

    pip install ObjectNat

----

Конфигурация
------------

Настройте вывод логов и прогресс-бары через модуль конфигурации::

    from objectnat import config

    config.change_logger_lvl("INFO")   # отключить отладочные логи
    config.set_enable_tqdm(False)      # отключить прогресс-бары tqdm

----

Контакты
--------

- `НЦКР <https://actcognitive.org/>`_ — Национальный центр когнитивных исследований
- `ИДУ <https://idu.itmo.ru/>`_ — Институт дизайна и урбанистики
- `Наталья Чичкова <https://t.me/nancy_nat>`_ — менеджер проекта
- `Данила Олейников (Donny) <https://t.me/ddonny_dd>`_ — ведущий инженер-разработчик

----

Публикации
----------

Скоро будут опубликованы.
