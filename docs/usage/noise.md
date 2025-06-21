# **Noise Simulation**

Noise Simulation models how sound propagates from one or more source points, taking into account **obstacles**, **vegetation**, and environmental conditions.  
The outputs are noise exposure maps useful for urban planning and environmental impact assessments.

---

The module provides two methods:

## **Full Wave-Based Simulation**

Performs detailed noise modeling using full wave-based calculations.

::: objectnat.simulate_noise
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/noise_simulation_1point.png" alt="noise_simulation_1point">


---

## **Simplified Noise Frame**

Creates a simplified noise exposure map using only geometric visibility and sound decay, without full wave simulation. Ideal for quick results.

::: objectnat.calculate_simplified_noise_frame
    options:
        show_root_heading: true
        heading_level: None

<img src="https://raw.githubusercontent.com/DDonnyy/ObjectNat/assets/noise_frame.png" alt="noise_frame">

---

## More Details

For comprehensive documentation and advanced configuration options, see the project Wiki:  
[Noise Simulation on GitHub](https://github.com/DDonnyy/ObjectNat/wiki/Noise-simulation)
