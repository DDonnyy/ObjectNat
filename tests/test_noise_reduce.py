"""Unit tests for the pure noise-attenuation helpers in ``noise_reduce``.

No geometry, graph, or data fixtures — just the interpolation table and the decay
formulas. Expected values are derived by hand from ``noise_init_data.air_resist_ratio``
(index = frequency Hz, columns = temperature °C).
"""

import pytest

from objectnat.methods.noise.noise_reduce import (
    dist_to_target_db,
    get_air_resist_ratio,
    green_noise_reduce_db,
)


# --------------------------------------------------------------------------- #
# get_air_resist_ratio
# --------------------------------------------------------------------------- #
def test_air_resist_ratio_exact_lookup():
    # Both temperature and frequency are tabulated exactly.
    assert get_air_resist_ratio(20, 2000) == pytest.approx(0.0096)
    assert get_air_resist_ratio(0, 8000) == pytest.approx(0.156)


def test_air_resist_ratio_interpolate_temperature_only():
    # freq = 500 exact; temp = 15 halfway between 10 (0.002) and 20 (0.0028).
    assert get_air_resist_ratio(15, 500) == pytest.approx(0.0024)


def test_air_resist_ratio_interpolate_frequency_only():
    # temp = 20 exact; freq = 750 halfway between 500 (0.0028) and 1000 (0.0052).
    assert get_air_resist_ratio(20, 750) == pytest.approx(0.004)


def test_air_resist_ratio_bilinear_interpolation():
    # temp = 15 (between 10/20), freq = 750 (between 500/1000): full bilinear blend.
    assert get_air_resist_ratio(15, 750) == pytest.approx(0.003475)


def test_air_resist_ratio_clamps_out_of_range():
    # temp below min (0) and freq above max (8000) clamp to the table corner loc[8000, 0].
    assert get_air_resist_ratio(-50, 100_000) == pytest.approx(0.156)


def test_air_resist_ratio_check_flag_still_returns_value():
    # check_temp_freq only emits warnings; the clamped value is still returned.
    assert get_air_resist_ratio(40, 63, check_temp_freq=True) == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# dist_to_target_db
# --------------------------------------------------------------------------- #
def test_dist_to_target_db_zero_air_resistance():
    # At 63 Hz the air resistance coefficient is 0, so decay is pure geometric
    # spreading: 20*log10(r) = init - target  ->  r = 10**((90-40)/20) = 10**2.5.
    dist = dist_to_target_db(90, 40, geometric_mean_freq_hz=63, air_temperature=20)
    assert dist == pytest.approx(10**2.5, rel=1e-4)


def test_dist_to_target_db_monotonic_in_target():
    # A higher target (closer to the source level) is reached over a shorter distance.
    far = dist_to_target_db(90, 40, geometric_mean_freq_hz=2000, air_temperature=20)
    near = dist_to_target_db(90, 60, geometric_mean_freq_hz=2000, air_temperature=20)
    assert far > near > 0


def test_dist_to_target_db_return_desc():
    desc = dist_to_target_db(90, 40, geometric_mean_freq_hz=2000, air_temperature=20, return_desc=True)
    assert isinstance(desc, str)
    assert "decays to 40 dB" in desc


# --------------------------------------------------------------------------- #
# green_noise_reduce_db
# --------------------------------------------------------------------------- #
def test_green_noise_reduce_db_known_value():
    # 0.08 * 10 * (8000**(1/3) / 8) = 0.08 * 10 * (20/8) = 2.0
    assert green_noise_reduce_db(8000, 10) == pytest.approx(2.0)


def test_green_noise_reduce_db_zero_thickness():
    assert green_noise_reduce_db(2000, 0) == pytest.approx(0.0)


def test_green_noise_reduce_db_monotonic_in_thickness():
    assert green_noise_reduce_db(2000, 20) > green_noise_reduce_db(2000, 5)
