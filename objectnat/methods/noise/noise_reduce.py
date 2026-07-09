import numpy as np
from scipy.optimize import fsolve

from objectnat import config

from .noise_init_data import air_resist_ratio

logger = config.logger


def get_air_resist_ratio(temp, freq, check_temp_freq=False):
    if check_temp_freq:
        if temp > max(air_resist_ratio.columns) or temp < min(air_resist_ratio.columns):
            logger.warning(
                f"The specified temperature of {temp}°C is outside the tabulated data range. "
                f"The air resistance coefficient for these values may be inaccurate. "
                f"Recommended temperature range: {min(air_resist_ratio.columns)}°C "
                f"to {max(air_resist_ratio.columns)}°C."
            )

        if freq > max(air_resist_ratio.index) or freq < min(air_resist_ratio.index):
            logger.warning(
                f"The specified geometric mean frequency of {freq} Hz is outside the tabulated data range."
                f" The air resistance coefficient for these values may be inaccurate."
                f" Recommended frequency range: {min(air_resist_ratio.index)} Hz to {max(air_resist_ratio.index)} Hz."
            )

    def get_nearest_values(array, value):
        sorted_array = sorted(array)
        if value in sorted_array:
            return [value]
        if value > max(sorted_array):
            return [sorted_array[-1]]
        if value < min(sorted_array):
            return [sorted_array[0]]

        for i, val in enumerate(sorted_array):
            if value < val:
                return sorted_array[max(i - 1, 0)], sorted_array[i]
        return sorted_array[-2], sorted_array[-1]

    nearest_temp = get_nearest_values(air_resist_ratio.columns, temp)
    nearest_freq = get_nearest_values(air_resist_ratio.index, freq)

    if len(nearest_temp) == 1 and len(nearest_freq) == 1:
        return air_resist_ratio.loc[nearest_freq[0], nearest_temp[0]]

    if len(nearest_temp) == 2 and len(nearest_freq) == 2:
        freq1, freq2 = nearest_freq
        temp1, temp2 = nearest_temp

        coef_temp1_freq1 = air_resist_ratio.loc[freq1, temp1]
        coef_temp1_freq2 = air_resist_ratio.loc[freq2, temp1]
        coef_temp2_freq1 = air_resist_ratio.loc[freq1, temp2]
        coef_temp2_freq2 = air_resist_ratio.loc[freq2, temp2]

        weight_temp1 = (temp2 - temp) / (temp2 - temp1)
        weight_temp2 = (temp - temp1) / (temp2 - temp1)
        weight_freq1 = (freq2 - freq) / (freq2 - freq1)
        weight_freq2 = (freq - freq1) / (freq2 - freq1)

        coef_freq1 = coef_temp1_freq1 * weight_temp1 + coef_temp2_freq1 * weight_temp2
        coef_freq2 = coef_temp1_freq2 * weight_temp1 + coef_temp2_freq2 * weight_temp2

        final_coef = coef_freq1 * weight_freq1 + coef_freq2 * weight_freq2

        return final_coef

    if len(nearest_temp) == 2 and len(nearest_freq) == 1:
        temp1, temp2 = nearest_temp
        freq1 = nearest_freq[0]

        coef_temp1 = air_resist_ratio.loc[freq1, temp1]
        coef_temp2 = air_resist_ratio.loc[freq1, temp2]

        weight_temp1 = (temp2 - temp) / (temp2 - temp1)
        weight_temp2 = (temp - temp1) / (temp2 - temp1)

        return coef_temp1 * weight_temp1 + coef_temp2 * weight_temp2

    if len(nearest_temp) == 1 and len(nearest_freq) == 2:
        temp1 = nearest_temp[0]
        freq1, freq2 = nearest_freq

        coef_freq1 = air_resist_ratio.loc[freq1, temp1]
        coef_freq2 = air_resist_ratio.loc[freq2, temp1]

        weight_freq1 = (freq2 - freq) / (freq2 - freq1)
        weight_freq2 = (freq - freq1) / (freq2 - freq1)

        return coef_freq1 * weight_freq1 + coef_freq2 * weight_freq2


def dist_to_target_db(
    init_noise_db, target_noise_db, geometric_mean_freq_hz, air_temperature, return_desc=False, check_temp_freq=False
) -> float | str:
    """
    Calculates the distance required for a sound wave to decay from an initial noise level to a target noise level,
    based on the geometric mean frequency of the sound and the air temperature. Optionally, can return a description
    of the sound propagation behavior.

    Args:
        init_noise_db (float): The initial noise level of the source in decibels (dB). This is the starting sound
            intensity.
        target_noise_db (float): The target noise level in decibels (dB), representing the level to which the sound
            decays over distance.
        geometric_mean_freq_hz (float): The geometric mean frequency of the sound (in Hz). This frequency influences
            the attenuation of sound over distance. Higher frequencies decay faster than lower ones.
        air_temperature (float): The temperature of the air in degrees Celsius. This influences the air's resistance
            to sound propagation.
        return_desc (bool, optional): If set to `True`, the function will return a description of the sound decay
            process instead of the calculated distance.
        check_temp_freq (bool, optional): If `True`, the function will check whether the temperature and frequency
            are within valid ranges.

    Returns:
        float or str: If `return_desc` is `False`, the function returns the distance (in meters) over which the sound
        decays from `init_noise_db` to `target_noise_db`. If `return_desc` is `True`, a descriptive string is returned
        explaining the calculation and the conditions.
    """

    def equation(r):
        return l - l_ist + 20 * np.log10(r) + k * r

    l_ist = init_noise_db
    l = target_noise_db
    k = get_air_resist_ratio(air_temperature, geometric_mean_freq_hz, check_temp_freq)
    initial_guess = 1
    r_solution = fsolve(equation, initial_guess)
    if return_desc:
        string = (
            f"Noise level of {init_noise_db} dB "
            f"with a geometric mean frequency of {geometric_mean_freq_hz} Hz "
            f"at an air temperature of {air_temperature}°C decays to {target_noise_db} dB "
            f"over a distance of {r_solution[0]} meters. Air resistance coefficient: {k}."
        )
        return string
    return r_solution[0]


def green_noise_reduce_db(geometric_mean_freq_hz, r_tree) -> float:
    """
    Calculates the amount of noise reduction (in dB) provided by vegetation of a given thickness at a specified
    geometric mean frequency. The function models the reduction based on the interaction of the sound with trees or
    vegetation.

    Args:
        geometric_mean_freq_hz (float): The geometric mean frequency of the sound (in Hz).
        r_tree (float): The thickness or density of the vegetation (in meters).

    Returns:
        float: The noise reduction (in dB) achieved by the vegetation. This value indicates how much quieter the sound
        will be after passing through or interacting with the vegetation of the specified thickness.
    """
    return round(0.08 * r_tree * ((geometric_mean_freq_hz ** (1 / 3)) / 8), 1)
