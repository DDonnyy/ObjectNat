import numpy as np


def min_max_normalization(data, new_min=0, new_max=1):
    """
    Min-max normalization for a given array of data.

    Parameters
    ----------
    data: numpy.ndarray
        Input data to be normalized.
    new_min: float, optional
        New minimum value for normalization. Defaults to 0.
    new_max: float, optional
        New maximum value for normalization. Defaults to 1.

    Returns
    -------
    numpy.ndarray
        Normalized data.

    Examples
    --------
    >>> import numpy as np
    >>> data = np.array([1, 2, 3, 4, 5])
    >>> normalized_data = min_max_normalization(data, new_min=0, new_max=1)
    """

    min_value = np.min(data)
    max_value = np.max(data)
    normalized_data = (data - min_value) / (max_value - min_value) * (new_max - new_min) + new_min
    return normalized_data
