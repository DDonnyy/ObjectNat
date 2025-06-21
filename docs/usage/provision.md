# **Service Provision Analysis**

This module evaluates how well services (e.g., schools, clinics) cover residential buildings based on their **capacity** and **accessibility**.  
It models demand-supply relationships and provides tools to analyze and adjust service coverage.

---

## **Evaluate initial provision**

Calculates provision scores between population points and service facilities considering:

- Distance or time threshold,
- Facility capacity,
- Demand distribution.

::: objectnat.get_service_provision
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_initial.png" alt="service_provision_initial">

---

## **Recalculate provision**

Allows you to recalculate provision results with new accessibility thresholds without recomputing the full OD-matrix.

::: objectnat.recalculate_links
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_recalculated.png" alt="service_provision_recalculated">

---

## **Clip to analysis area**

Restricts provision output to a given geographic boundary (e.g., administrative area).

::: objectnat.clip_provision
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/service_provision_clipped.png" alt="service_provision_clipped">

---

## Graph Preparation via IduEdu

For best performance and reproducibility, it is recommended to build or download intermodal graphs using the [`IduEdu`](https://pypi.org/project/iduedu/) library.

Here is a minimal example:

```python
# Install required packages (uncomment if needed)
# !pip install iduedu

from iduedu import get_boundary, get_intermodal_graph

# Load boundary and graph for a specific region using OSM ID 1114252.
poly = get_boundary(osm_id=1114252)
G_intermodal = get_intermodal_graph(polygon=poly, clip_by_bounds=True)
```