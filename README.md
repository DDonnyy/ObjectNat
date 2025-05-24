# ObjectNat

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://img.shields.io/pypi/v/objectnat.svg)](https://pypi.org/project/objectnat/)
[![CI](https://github.com/DDonnyy/ObjectNat/actions/workflows/ci_pipeline.yml/badge.svg)](https://github.com/DDonnyy/ObjecNat/actions/workflows/ci_pipeline.yml)
[![codecov](https://codecov.io/gh/DDonnyy/ObjectNat/graph/badge.svg?token=K6JFSJ02GU)](https://codecov.io/gh/DDonnyy/ObjectNat)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://opensource.org/licenses/MIT)

- [РИДМИ (Russian)](README_ru.md)
<p align="center">
<img src="https://github.com/user-attachments/assets/ed0f226e-1728-4659-9e21-b4d499e703cd" alt="logo" width="400">
</p>

#### **ObjectNat** is an open-source library created for geospatial analysis created by **IDU team**

## Features and how to use

1. **[Isochrones and Transport Accessibility](./examples/isochrone_generator.ipynb)** — Isochrones represent areas reachable from a starting point within a given time limit along a transport network. This function enables analysis of transport accessibility using pedestrian, automobile, public transport graphs, or their combination.

   The library offers multiple isochrone generation methods:
   - **Baseline isochrones**: show a single area reachable within a specified time.
   - **Stepped isochrones**: show accessibility ranges divided into time intervals (e.g., 5, 10, 15 minutes).

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

2. **[Coverage Zones](./examples/coverage_zones.ipynb)** — Function for generating **coverage zones** from a set of source points using a transport network. It calculates the area each point can reach based on **travel time** or **distance**, then builds polygons via **Voronoi diagrams** and clips them to a custom boundary if provided.

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_time_10min.png" alt="coverage_zones_time_10min" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_distance_600m.png" alt="coverage_zones_distance_600m" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_radius_distance_800m.png" alt="coverage_zones_distance_radius_voronoi" width="350">
   </p>
   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_separate.png" alt="stepped_coverage_zones_separate" width="350">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_voronoi.png" alt="stepped_coverage_zones_voronoi" width="350">
   </p>
 
3. **[Service Provision Analysis](./examples/calculate_provision.ipynb)** — Function for evaluating the provision of residential buildings and their population with services (e.g., schools, clinics)
    that have limited **capacity** and a defined **accessibility threshold** (in minutes or distance). The function models **demand-supply balance**, estimating how well services meet the needs of nearby buildings within the allowed time.

   The library also supports:
   - **Recalculation** of existing provision results using a new time threshold.
   - **Clipping** of provision results to a custom analysis area (e.g., administrative boundaries).

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_initial.png" alt="service_provision_initial" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_recalculated.png" alt="service_provision_recalculated" width="300">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_clipped.png" alt="service_provision_clipped" width="300">
   </p>

4. **[Visibility Analysis](./examples/visibility_analysis.ipynb)** — Function for estimating visibility from a given point or multiple points to nearby buildings within a certain distance.
   This can be used to assess visual accessibility in urban environments.
   The library also includes a **catchment area calculator** for large-scale visibility analysis based on a dense grid of observer points (recommended: ~1000 points spaced 10–20 meters apart).
   Points can be generated using a road network and distributed along edges.

   The module includes:
   - A **fast approximate method** for large datasets.
   - A **accurate method** for detailed local analysis.

   <p align="center">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/visibility_comparison_methods.png" alt="visibility_comparison_methods" height="250">
     <img src="https://github.com/user-attachments/assets/b5b0d4b3-a02f-4ade-8772-475703cd6435" alt="visibility-catchment-area" height="250">
   </p>

5. **[Noise Simulation](./examples/noise_simulation.ipynb)** — Simulates noise propagation from a set of source points, taking into account **obstacles**, **vegetation**, and **environmental factors**.

   🔗 **[See detailed explanation in the Wiki](https://github.com/DDonnyy/ObjectNat/wiki/Noise-simulation)**

   <p align="center">
      <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/noise_simulation_1point.png" alt="noise_simulation_1point" width="400">
   </p>


6. **[Point Clusterization](./examples/point_clusterization.ipynb)** — Function to generate **cluster polygons** from a set of input points based on:
   - Minimum **distance** between points.
   - Minimum **number of points** per cluster.

   Additionally, the function can calculate the **relative ratio** of different service types within each cluster, enabling spatial analysis of service composition.

   <p align="center">
     <img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/building_clusters.png" alt="building_clusters" width="400">
   </p>

## City graphs

To ensure optimal performance of ObjectNat's geospatial analysis functions, it's recommended to utilize urban graphs sourced from the [IduEdu](https://github.com/DDonnyy/IduEdu) library.
**IduEdu** is an open-source Python library designed for the creation and manipulation of complex city networks derived from OpenStreetMap data. 

**IduEdu** can be installed with ``pip``:
```
pip install IduEdu
```
## Installation

**ObjectNat** can be installed with ``pip``:

```
pip install ObjectNat
```
### Configuration changes

```python
from objectnat import config

config.change_logger_lvl('INFO')  # To mute all debug msgs
config.set_enable_tqdm(False)  # To mute all tqdm's progress bars
```
## Contacts

- [NCCR](https://actcognitive.org/) - National Center for Cognitive Research
- [IDU](https://idu.itmo.ru/) - Institute of Design and Urban Studies
- [Natalya Chichkova](https://t.me/nancy_nat) - project manager
- [Danila Oleynikov (Donny)](https://t.me/ddonny_dd) - lead software engineer

## Publications
