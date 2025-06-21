# **Isochrones and Transport Accessibility**

Isochrones represent areas reachable from a starting point within a given time limit along a transport network.
This function enables analysis of transport accessibility using pedestrian, automobile, public transport graphs, or their combination.

---

The library offers multiple isochrone generation methods:

## **Baseline isochrones**
Show a single area reachable within a specified time.

::: objectnat.get_accessibility_isochrones
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_ways_15_min.png" alt="isochrone_ways_15_min">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_radius_15_min.png" alt="isochrone_radius_15_min">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/isochrone_3points_radius_8_min.png" alt="isochrone_3points_radius_8_min">

---

## **Stepped isochrones**
Show accessibility ranges divided into time intervals (e.g., 5, 10, 15 minutes).

::: objectnat.get_accessibility_isochrone_stepped
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_ways_15_min.png" alt="stepped_isochrone_ways_15_min">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_radius_15_min.png" alt="stepped_isochrone_radius_15_min">
<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/stepped_isochrone_separate_15_min.png" alt="stepped_isochrone_separate_15_min">

---