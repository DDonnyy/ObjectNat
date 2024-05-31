# ObjectNat - Meta Library

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

<p align="center">
<img src="https://i.ibb.co/FWtHNQv/logo.png" alt="logo" width="400">
</p>

#### **ObjectNat** is an open-source library created for geospatial analysis created by **IDU team**

## ObjectNat Components

- [dongraphio](https://github.com/DDonnyy/dongraphio) : `dongraphio` provides graph functions
- [provisio](https://github.com/DDonnyy/provisio) : `provisio` provides main provisio fuctions
- [population-restorator](https://github.com/kanootoko/population-restorator) : `restorator` provides city resettlement

## Features and how to use

1. **[City graph from OSM](./examples/graph_generator.ipynb)** - Function to assemble a road, pedestrian, and public
   transport graph from OpenStreetMap (OSM) and creating Intermodal graph.

   <img src="https://i.ibb.co/VpsPgL1/Graph-generator-v1.webp" alt="Graph-generator-v1" height="250">

2. **[Adjacency matrix](./examples/calculate_adjacency_matrix.ipynb)** - Calculate adjacency matrix based on the provided
   graph and edge weight type (time or distance). The intermodal graph can be obtained using the previous example.
3. **[Isochrones,transport accessibility](./examples/isochrone_generator.ipynb)** - Function for generating isochrones to
   analyze transportation accessibility from specified starting coordinates. Isochrones can be constructed based on
   pedestrian, automobile, or public transport graphs, or a combination thereof.

   <img src="https://i.ibb.co/QM0tmZ2/isochrones-from-2-points.webp" alt="isochrones-from-2-points" height="250">

4. **[Population restoration](./examples/restore_population.ipynb)** - Function for resettling population into the provided
   layer of residential buildings. This function distributes people among dwellings based on the total city population
   and the living area of each house.
5. **[Service provision](./examples/calculate_provision.ipynb)** - Function for calculating the provision of residential
   buildings and population with services. In case of missing data, this function utilizes previously described
   functionality to retrieve the necessary information.

   <img src="https://i.ibb.co/CW7Xj5F/Burger-Provision5min.webp" alt="Burger-Provision5min" height="250">
   
6. **[Visibility analysis](./examples/visibility_analysis.ipynb)** - Function to get a quick estimate of visibility from a
   given point(s) to buildings within a given distance. Also, there is a visibility catchment area calculator for a
   large
   urban area. This function is designed to work with at least 1000 points spaced 10-20 meters apart for optimal
   results. Points can be generated using a road graph and random point distribution along edges.

   <img src="https://i.ibb.co/LxcGTfN/visibility-from-point.webp" alt="visibility-from-point" height="250"> 

   <img src="https://i.ibb.co/zNRzXc5/visibility-catchment-area.webp" alt="visibility-catchment-area" height="250">

7. **[Point clusterization](./examples/point_clusterization.ipynb)** - Function to generate cluster polygons for given
   points based on a specified minimum distance and minimum points per cluster. Optionally, calculate the relative ratio
   between types of services within the clusters.

   <img src="https://i.ibb.co/fF3c4YC/service-clusterization.webp" alt="service-clusterization" height="250">

## Installation

**ObjectNat** can be installed with ``pip``:

```
pip install ObjectNat
```

## Contacts

- [NCCR](https://actcognitive.org/) - National
  Center for Cognitive Research
- [IDU](https://idu.itmo.ru/) - Institute of
  Design and Urban Studies
- [Natalya Chichkova](https://t.me/nancy_nat) - project manager
- [Danila Oleynikov (Donny)](https://t.me/ddonny_dd) - lead software engineer

## Publications
