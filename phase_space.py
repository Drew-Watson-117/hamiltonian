"""
This file contains a class that represents a phase space
"""
from typing import Optional
from sympy import Symbol, Function, Eq

class PhaseSpace:
    """
    Class representing phase space
    """
    
    def __init__(self, canonical_coordinates: dict[Symbol, Symbol], independent_variable: Symbol, hamiltonian: Optional[Function] = None):
        """
        Construct a phase space given a set of canonical coordinates

        Parameters:
            canonical_coordinates (dict[Symbol, Symbol]):
                A dictionary containing the canonical coordinates
                with the positions as the keys and the momenta as the values
            independent_variable (Symbol):
                The independent variable for the space (time for classical mechanics)
            hamiltonian (Function):
                The hamiltonian function for the phase space, if there is one
        """
        self.hamiltonian = hamiltonian
        self.independent_variable = independent_variable
        self.positions = canonical_coordinates.keys()
        self.momenta = canonical_coordinates.values()
        self.dimension = len(canonical_coordinates)

        self._canonical_coordinates = canonical_coordinates
        self._momenta_to_positions = {}
        for pos, mom in self._canonical_coordinates.items():
            self._momenta_to_positions[mom] = pos

    def set_hamiltonian(self, hamiltonian: Function):
        """
        Set the Hamiltonian function for the phase space
        """
        self.hamiltonian = hamiltonian

    def get_momentum(self, position: Symbol):
        """
        Get the conjugate momentum of a given position
        """
        if position not in self._canonical_coordinates:
            raise ValueError(f"Position {position} not in the phase space")
        return self._canonical_coordinates[position]

    def get_position(self, momentum: Symbol):
        """
        Get the conjugate position of a given momentum
        """
        if momentum not in self._momenta_to_positions:
            raise ValueError(f"Position {momentum} not in the phase space")
        return self._momenta_to_positions[momentum]

    def hamiltons_equation(self, *coordinates: Symbol) -> list[Eq]:
        """
        Given a list of coordinates, return Hamilton's Equation 
        for that position or momentum
        """
        H = self._substitute_functions_in_hamiltonian()
        t = self.independent_variable
        eqns = []
        for pos_or_mom in coordinates:
            if pos_or_mom in self.positions:
                x = self._convert_to_function(pos_or_mom)
                p = self._convert_to_function(self.get_momentum(pos_or_mom))
                eqns.append(Eq(x.diff(t),H.diff(p)))
            elif pos_or_mom in self.momenta:
                p = self._convert_to_function(pos_or_mom)
                x = self._convert_to_function(self.get_position(pos_or_mom))
                eqns.append(Eq(p.diff(t), -H.diff(x)))
            else:
                raise ValueError(f"Coordinate {pos_or_mom} not found in the space")
        return eqns

    def hamiltons_equations(self) -> list[Eq]:
        """
        Calculate Hamilton's Equations for the space

        Returns:
            Hamilton's Equatiosn (list[Eq])
        """
        eqns: list[Eq] = []
        H = self._substitute_functions_in_hamiltonian()
        t = self.independent_variable
        for pos, mom in self._canonical_coordinates.items():
            x = self._convert_to_function(pos)
            p = self._convert_to_function(mom)
            pos_eqn = Eq(
                x.diff(t),
                H.diff(p)
            )
            mom_eqn = Eq(
                p.diff(t),
                -H.diff(x)
            )
            eqns.extend([pos_eqn, mom_eqn])
        return eqns


    def _convert_to_function(self, sym: Symbol) -> Function:
        """
        Convert the given Symbol to a function of the independent
        variable
        """
        return Function(f"{sym}")(self.independent_variable)

    def _substitute_functions_in_hamiltonian(self) -> Function:
        """
        Substitute the symbols for dependent variables
        to functions in the Hamiltonian
        """
        H = self.hamiltonian
        for pos, mom in self._canonical_coordinates.items():
            x = self._convert_to_function(pos)
            p = self._convert_to_function(mom)
            # Substitute the position variable for the position function
            H = H.subs(pos, x)
            # Substitute the momentum variable for the momentum function
            H = H.subs(mom, p)
        return H

    def __repr__(self):
        pairs = [f"({pos}, {mom})" for pos, mom in self._canonical_coordinates.items()]
        return " ".join(pairs)
