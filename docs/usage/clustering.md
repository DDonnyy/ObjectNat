# **Point Clusterization**

Clusterization groups nearby points into polygons based on spatial proximity and minimum cluster size.  
It can be used for identifying dense urban areas, service hubs, or catchment zones.

---

## **Cluster Generation**

Clusters are generated using spatial rules:

- **Minimum distance** between points to be included in the same cluster.
- **Minimum number of points** required to form a valid cluster.

::: objectnat.get_clusters_polygon
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/building_clusters.png" alt="building_clusters">

---

