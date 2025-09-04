"""
This file contains the Ostrogradsky algorithm
for computing Hamiltonians
"""

from sympy import Symbol, Function, Matrix, Eq, solve, simplify, expand
from .coordinate_system import CoordinateSystem
from .configuration_space import ConfigurationSpace, DependentVariable
from .phase_space import PhaseSpace

def _hessian_determinant(f: Function, coordinates: list[Symbol]):
    """
    Given a function and coordinate system, compute the determinant
    of the Hessian matrix
    """
    hessian_matrix = Matrix([
        [f.diff(x).diff(y) for x in coordinates] 
        for y in coordinates])
    return hessian_matrix.det()

def hamiltonian(L: Function, coordinate_system: CoordinateSystem):
    """
    Given a Lagrangian L, a system of coordinates, and a list of
    the highest order derivative present in the Lagrangian for each coordinate,
    calculate the Hamiltonian using Ostrogradsky's Formalism

    Parameters:
        L (Function): 
            The Lagrangian
        coordinate_system (CoordinateSystem): 
            The coordinate system defining coordinates in the Lagrangian
    Returns:
        phase space with the canonical coordinates and the Hamiltonian (PhaseSpace)
    """
    

    space = ConfigurationSpace(L, coordinate_system)
    L = space.L

    # Define "X" as the dependent variables (and their derivatives)
    X = space.dependent_variables
    highest_order_derivatives = [x[-1] for x in X.values()]
    # Check that the Hessian is non-degenerate
    if _hessian_determinant(L, highest_order_derivatives) == 0:
        raise ValueError("The Hessian of L must be non-degenerate "
                         "(determinant of the Hessian matrix must be non-zero)")

    # Find A for each dependent variable
    A_s: dict[DependentVariable, Function] = {}
    canonical_coords: dict[Symbol, Symbol] = {}
    for x in X.values(): # labels which dependent variable is considered
        p_equations: list[Eq] = []
        # N = len(x)
        N = len(x) - 1
        # for i in range(N-1):
        for i in range(N):
            P = Symbol("P_{" + str(x[i]) + "}")
            canonical_coords[x[i]] = P
            rhs = sum((-1)**(j-i-1) * space.total_derivative(L.diff(x[j]), j-i-1) for j in range(i+1, N+1))
            p_equation = Eq(P, rhs)
            p_equations.append(p_equation)
        sln = solve(p_equations, x[-1])
        A_s = A_s | sln
    phase_space = PhaseSpace(
        canonical_coords, 
        independent_variable=coordinate_system.independent_variable
    )
    # Compute Hamiltonian
    H = 0
    for x in X.values():
        term = 0
        for i in range(N-1):
            x_current = x[i]
            x_dot = x[i+1]
            term += phase_space.get_momentum(x_current) * x_dot
        p_last = phase_space.get_momentum(x[N-1])
        # Add "A" term
        H += term + p_last * A_s[x[N]]
    # Substitute velocities for momenta in the Lagrangian
    L_p = L
    for vel, mom in A_s.items():
        L_p = L_p.subs(vel, mom)
    # Subtract the Lagrangian to get the Hamiltonian
    H -= L_p
    phase_space.set_hamiltonian(expand(simplify(H)))
    # phase_space.set_hamiltonian(H)
    return phase_space