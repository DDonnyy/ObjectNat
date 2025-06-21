# **Visibility Analysis**

Visibility analysis estimates which buildings or areas are visible from a given observer point (or set of points)
within a specific distance. This is useful in assessing **visual accessibility**, urban form,
and perceptual exposure in public space.

---

The module supports multiple modes of analysis:

## **Accurate Method**

Computes visibility using fine-grained raster-based methods. More accurate for local areas, but slower.

::: objectnat.get_visibility_accurate
    options:
        show_root_heading: true
        heading_level: None

---

## **Fast Approximate Method**

Optimized for large datasets or large areas. Uses geometry simplifications and vector-based visibility.

::: objectnat.get_visibility
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/visibility_comparison_methods.png" alt="visibility_comparison_methods">

---

## **Catchment Visibility from Multiple Points**

Performs visibility analysis for a dense grid of observer points.  
Used to generate **catchment areas** of visibility (e.g., “where can this building be seen from?”).

::: objectnat.get_visibilities_from_points
    options:
        show_root_heading: true
        heading_level: None

---

_The image below shows an example of using visibility polygons to calculate "visibility pools" - areas in an urban 
environment that are most visible from different locations._

<img src="https://github.com/user-attachments/assets/b5b0d4b3-a02f-4ade-8772-475703cd6435" alt="visibility-catchment-area">

---
