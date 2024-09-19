# ObjectNat - Meta Library

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://img.shields.io/pypi/v/objectnat.svg)](https://pypi.org/project/objectnat/)

- [РИДМИ (Russian)](README_ru.md)
<p align="center">
<img src="https://github.com/user-attachments/assets/d3878cce-8eba-4f96-8458-9a798d436120" alt="logo" width="400">
</p>

#### **ObjectNat** is an open-source library created for geospatial analysis created by **IDU team**

## ObjectNat Components

- [IduEdu](https://github.com/DDonnyy/IduEdu) : `IduEdu` provides graph functions
- [population-restorator](https://github.com/kanootoko/population-restorator) : `restorator` provides city resettlement

## Features and how to use

1. **[City graph from OSM (IduEdu)](https://github.com/DDonnyy/IduEdu/blob/main/examples/get_any_graph.ipynb)** - Functions to assemble a road, pedestrian,
   and public transport graph from OpenStreetMap (OSM) and creating Intermodal graph. 

   <img src="https://github.com/user-attachments/assets/8dc98da9-8462-415e-8cc8-bdfca788e206" alt="IntermodalGraph" height="250">

2. **[Adjacency matrix](./examples/calculate_adjacency_matrix.ipynb)** - Calculate adjacency matrix based on the provided
   graph and edge weight type (time or distance). The intermodal graph can be obtained using the previous example.

3. **[Isochrones,transport accessibility](./examples/isochrone_generator.ipynb)** - Function for generating isochrones to
   analyze transportation accessibility from specified starting coordinates. Isochrones can be constructed based on
   pedestrian, automobile, or public transport graphs, or a combination thereof.

   <img src="https://github.com/user-attachments/assets/37f308a5-db56-497d-b080-4edef3584fe5" alt="isochrones" height="250">

4. **[Population restoration](./examples/restore_population.ipynb)** - Function for resettling population into the provided
   layer of residential buildings. This function distributes people among dwellings based on the total city population
   and the living area of each house.
5. **[Service provision](./examples/calculate_provision.ipynb)** - Function for calculating the provision of residential
   buildings and population with services. 

   <img src="https://github.com/user-attachments/assets/5f2b3c55-9a02-4d70-80f4-503b77023eda" alt="ProvisionSchools" height="250">

6. **[Visibility analysis](./examples/visibility_analysis.ipynb)** - Function to get a quick estimate of visibility from a
   given point(s) to buildings within a given distance. Also, there is a visibility catchment area calculator for a
   large
   urban area. This function is designed to work with at least 1000 points spaced 10-20 meters apart for optimal
   results. Points can be generated using a road graph and random point distribution along edges.

   <img src="https://github.com/user-attachments/assets/2927ac86-01e8-4b0e-9ea8-72ad81c13cf5" alt="visibility-from-point" height="250"> 

   <img src="https://github.com/user-attachments/assets/b5b0d4b3-a02f-4ade-8772-475703cd6435" alt="visibility-catchment-area" height="250">
   
7. **[Point clusterization](./examples/point_clusterization.ipynb)** - Function to generate cluster polygons for given
   points based on a specified minimum distance and minimum points per cluster. Optionally, calculate the relative ratio
   between types of services within the clusters.

   <img src="https://github.com/user-attachments/assets/2a9ad722-87d2-4954-9612-5ac3765aa824" alt="service-clusterization" height="250">
   
8. **[Living buildings from OSM](./examples/download_buildings_from_osm.ipynb)** - This function downloads building geometries from OpenStreetMap (OSM) for a specified territory and assigns attributes to each building. Specifically, it determines whether a building is residential (`is_living` attribute) and estimates the approximate number of inhabitants (`approximate_pop` attribute).
   
   <img src="https://github.com/user-attachments/assets/d60dcd85-1a2e-4342-aae4-561aeda18858" alt="Living buildings" height="250">
   
## Installation

**ObjectNat** can be installed with ``pip``:

```
pip install ObjectNat
```
### Configuration changes

```python
from objectnat import config

config.set_timeout(10)  # Timeout for overpass queries
config.change_logger_lvl('INFO')  # To mute all debug msgs
config.set_enable_tqdm(False)  # To mute all tqdm's progress bars
config.set_overpass_url('http://your.overpass-api.de/interpreter/URL')
```
## Contacts

- [NCCR](https://actcognitive.org/) - National
  Center for Cognitive Research
- [IDU](https://idu.itmo.ru/) - Institute of
  Design and Urban Studies
- [Natalya Chichkova](https://t.me/nancy_nat) - project manager
- [Danila Oleynikov (Donny)](https://t.me/ddonny_dd) - lead software engineer

## Publications
