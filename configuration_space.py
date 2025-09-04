from sympy import Symbol, Function
from sympy.solvers.deutils import ode_order

from .coordinate_system import CoordinateSystem

class DependentVariable(Symbol):
    """
    Class representing a dependent variable
    """

    @classmethod
    def from_function(cls, fcn: Function) -> "DependentVariable":
        """
        Create a dependent variable with the same name
        as the given function
        """
        return cls(fcn.name)

    def total_diff(self, n: int = 1) -> "DependentVariable":
        """
        Compute the n-th derivative of a dependent variable
        """
        if n < 0:
            raise ValueError("Cannot compute a negative-order derivative")
        elif n == 0:
            return self
        underscore_split = str(self).split("_")
        if underscore_split[-1].isnumeric():
            i = int(underscore_split.pop(-1))
            underscore_split.append(str(i+n))
        else:
            underscore_split.append(str(n))
        name = "_".join(underscore_split)
        return DependentVariable(name)

class ConfigurationSpace:
    """
    Class that represents the configuration space
    (coordinates, the lagrangian, and all derivatives of the Lagrangian)
    """
    def __init__(
            self,
            L: Function,
            coordinate_system: CoordinateSystem):
        """
        Create a CoordinatesAndDerivatives object
        """
        self.L = L
        self.coordinate_system = coordinate_system
        self.independent_variable: Symbol = coordinate_system.independent_variable
        self.dependent_variables: dict[str, list[DependentVariable]] = {var_name: [] for var_name in coordinate_system.dependent_variables.keys()}
        self.dimension = coordinate_system.dimension
        
        t = self.independent_variable
        for x_name, x_fcn in coordinate_system.dependent_variables.items():
            order = ode_order(L, x_fcn)
            x = DependentVariable.from_function(x_fcn)
            # Substitute the function for the coordinate
            for n in range(order, -1, -1):
                x_n = x.total_diff(n)
                x_n_fcn = x_fcn.diff(t, n)
                self.L = self.L.subs(x_n_fcn, x_n)
                self.dependent_variables[x_name].append(x_n)
            self.dependent_variables[x_name].reverse()
            

    def get_coordinate(self, name: str) -> DependentVariable:
        """
        Given the name of a coordinate, return the coordinate
        """
        if name in self.dependent_variables:
            return self.dependent_variables[name][0]
        else:
            raise ValueError(f"Coordinate {name} is not in the space")

    def get_coordinates(self) -> list[DependentVariable]:
        """
        Return a list of the coordinates in the space
        """
        coords = []
        for coords_and_derivatives in self.dependent_variables.values():
            coords.append(coords_and_derivatives[0])
        return coords
    
    def get_all_coordinates(self) -> list[DependentVariable]:
        """
        Return the list of all coordinates and their derivatives
        """
        # coords: list[JetVariable] = []
        coords: list[Symbol] = []
        for coord_list in self.dependent_variables.values():
            coords.extend(coord_list)
        return coords

    def total_derivative(self, f: Function, n: int = 1) -> Function:
        """
        Given a function f, compute the n-th total
        derivative of the function
        """
        if n < 0:
            raise ValueError("Cannot compute a negative order derivative")
        if n == 0:
            return f
        coords = self.get_all_coordinates()
        terms = [f.diff(x) * x.total_diff() for x in coords]
        result = sum(terms) + f.diff(self.independent_variable)
        return self.total_derivative(result, n-1)