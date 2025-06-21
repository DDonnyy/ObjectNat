# **Coverage Zones**

Coverage zones show areas that can be reached from each of multiple source points within a certain time or distance limit using a transport network.  
They are built by calculating reachability per point, generating Voronoi polygons, and optionally clipping them to a boundary.

---

The library supports several methods for generating coverage zones:

## **Coverage using transport graph**

Uses a full routing engine to determine reachable areas per point, then builds zones.

::: objectnat.get_graph_coverage
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_time_10min.png" alt="coverage_zones_time_10min">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_distance_600m.png" alt="coverage_zones_distance_600m">

---

## **Coverage using radius only**

Generates fixed-radius buffers per point without routing, clipped via Voronoi.

::: objectnat.get_radius_coverage
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/coverage_zones_radius_distance_800m.png" alt="coverage_zones_distance_radius_voronoi">

---

## **Stepped graph coverage**

Creates stepped zones (e.g., 5, 10, 15 minutes) using the full transport graph per point.

::: objectnat.get_stepped_graph_coverage
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_separate.png" alt="stepped_coverage_zones_separate">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_coverage_zones_voronoi.png" alt="stepped_coverage_zones_voronoi">

---

