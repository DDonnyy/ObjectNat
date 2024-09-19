class CapacityKeyError(KeyError):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "CapacityKeyError, {0} ".format(self.message)

        return (
            "Column 'capacity' was not found in provided 'services' GeoDataFrame. This attribute "
            "corresponds to the total capacity for each service."
        )


class CapacityValueError(ValueError):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "CapacityValueError, {0} ".format(self.message)

        return "Column 'capacity' in 'services' GeoDataFrame  has no valid value."


class DemandKeyError(KeyError):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "DemandKeyError, {0} ".format(self.message)

        return (
            "The column 'demand' was not found in the provided 'demanded_buildings' GeoDataFrame. "
            "This attribute corresponds to the number of demands for the selected service in each building."
        )


class DemandValueError(ValueError):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "DemandValueError, {0} ".format(self.message)
        return "Column 'demand' in 'demanded_buildings' GeoDataFrame  has no valid value."
