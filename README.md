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

---

## Features and How to Use

Each feature is accompanied by a Jupyter notebook example and full documentation.

1. **[Isochrones and Transport Accessibility](./examples/isochrone_generator.ipynb)**  
   Analyze areas reachable within a given time along a transport network.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/isochrones.html)

2. **[Coverage Zones](./examples/coverage_zones.ipynb)**  
   Build zones of reachability for each point using routing or simple radius.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/coverage.html)

3. **[Service Provision Analysis](./examples/calculate_provision.ipynb)**  
   Evaluate service availability and model demand-supply balance.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/provision.html)

4. **[Visibility Analysis](./examples/visibility_analysis.ipynb)**  
   Estimate visibility to nearby buildings from selected points.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/visibility.html)

5. **[Noise Simulation](./examples/noise_simulation.ipynb)**  
   Simulate noise propagation considering obstacles and environment.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/noise.html)  
   🔗 [Detailed theory in the Wiki](https://github.com/DDonnyy/ObjectNat/wiki/Noise-simulation)

6. **[Point Clusterization](./examples/point_clusterization.ipynb)**  
   Group nearby points into clusters and analyze service composition.  
   📄 [See documentation](https://iduclub.github.io/ObjectNat/latest/usage/clustering.html)

---

## City graphs

To ensure optimal performance of ObjectNat's geospatial analysis functions, it's recommended to utilize urban graphs sourced from the [IduEdu](https://github.com/DDonnyy/IduEdu) library.
**IduEdu** is an open-source Python library designed for the creation and manipulation of complex city networks derived from OpenStreetMap data. 

**IduEdu** can be installed with ``pip``:
```
pip install IduEdu
```
---

## Installation

**ObjectNat** can be installed with ``pip``:

```
pip install ObjectNat
```

---

### Configuration changes

```python
from objectnat import config

config.change_logger_lvl('INFO')  # To mute all debug msgs
config.set_enable_tqdm(False)  # To mute all tqdm's progress bars
```

---

## Contacts

- [NCCR](https://actcognitive.org/) - National Center for Cognitive Research
- [IDU](https://idu.itmo.ru/) - Institute of Design and Urban Studies
- [Natalya Chichkova](https://t.me/nancy_nat) - project manager
- [Danila Oleynikov (Donny)](https://t.me/ddonny_dd) - lead software engineer

---

## Publications

_Coming soon._