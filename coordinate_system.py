"""
The file contains the JetVariable and JetSpace
classes
"""

from sympy import symbols, Symbol, Function

class CoordinateSystem:
    """
    Defines a coordinate system
    """
    def __init__(
            self,
            independent_variable: str, 
            dependent_variables: list[str]):
        """
        Create a coordinate system
        """
        self.independent_variable: Symbol = symbols(independent_variable)
        # self.dependent_variables: dict[str, DependentVariable] = {var: DependentVariable(var) for var in dependent_variables}
        self.dependent_variables: dict[str, Function] = {var: Function(var)(self.independent_variable) for var in dependent_variables}
        self.dimension = len(self.dependent_variables)

    def get_coordinate(self, name: str) -> Function:
        """
        Given the name of a coordinate, return the coordinate
        """
        if name in self.dependent_variables:
            return self.dependent_variables[name]
        else:
            raise ValueError(f"Coordinate {name} is not in the space")

    def get_coordinates(self) -> list[Function]:
        """
        Return a list of the coordinates in the space
        """
        return list(self.dependent_variables.values())
