# `hamiltonian` Package

This document outlines the requirements for and usage of the `hamiltonian` package.

The purpose of the `hamiltonian` package is to conduct Hamiltonian mechanics given
a Lagrangian and a coordinate system.
This package computes Hamiltonians for Lagrangians that depend on higher order derivatives of the coordinates ($\ddot{x}, x^{(3)},...$).
This package uses `sympy` for symbolic manipulation.

## Requirements

This script requires the following:

- Python 3.11.0 or later
- `sympy` 1.11.1 or later

As long as you have a version of Python >= 3.11.0, you can install `sympy` using the `requirements.txt` file supplied with this code.
This can be done running the following command:

```bash
pip3 install -r requirements.txt
```

from **within this directory**.

## Usage

You can use the `hamiltonian` package by copying this directory into your desired code-base and importing it:

```python
import hamiltonian
```

## Defining a Coordinate System

The `hamiltonian` package includes a `CoordinateSystem` class.
This class is used to define the coordinate system used in your classical mechanics problem.
We supply the independent variable to the constructor as a string and the dependent variables to the constructor as a list of strings.
We can define 3 dimensional cartesian coordinates like this:

```python
import hamiltonian
cartesian_coordinates = hamiltonian.CoordinateSystem("t", ["x", "y", "z"])
```

We can get the independent and dependent variables this way:

```python
t = cartesian_coordinates.independent_variable
x, y, z = cartesian_coordinates.get_coordinates()
# Alternatively:
x = cartesian_coordinates.get_coordinate("x")
y = cartesian_coordinates.get_coordinate("y")
z = cartesian_coordinates.get_coordinate("z")
```

Here, `t` is a `sympy.Symbol` object and `x`, `y`, `z` are `sympy.Function` objects.

## Example 1: Computing the Hamiltonian for a Particle in a Uniform Gravitational Field

This example computes the Hamiltonian for a particle in a uniform gravitational field. The Lagrangian for the system is given by:

$$\mathcal{L} = \frac{1}{2}m(\dot{x}^2+ \dot{y}^2 + \dot{z}^2) - mgz$$

We can define this Lagrangian using `hamiltonian` and `sympy`  like-so:

```python
import hamiltonian
from sympy import symbols, Rational, init_printing(), Eq

# Define the coordinate system and constants
coordinate_system = hamiltonian.CoordinateSystem("t", ["x", "y", "z"])
t = coordinate_system.independent_variable
x, y, z = coordinate_system.get_coordinates()
m, g = symbols(["m", "g"])

# Define the Lagrangian
L = Rational(1,2) * m * (x.diff(t)**2 + y.diff(t)**2 + z.diff(t)**2) - m*g*z

# Print the Lagrangian
init_printing()
print(Eq("\mathcal{L}", L))
```

Note how derivatives of `x`, `y`, and `z` are obtained by calling `x.diff(t)`.

The Hamiltonian can be computed by calling the `hamiltonian.hamiltonian` function:

```python
phase_space = hamiltonian.hamiltonian(L, coordinate_system)
```

Note that the `hamiltonian.hamiltonian` function returns a `PhaseSpace` object. The `PhaseSpace` object contains the Hamiltonian:

```python
H = phase_space.hamiltonian
# Print the Hamiltonian
print(Eq("\mathcal{H}", H))
```

The phase space objects also supplies attributes and methods for manipulating positions and momenta:

```python
x, y, z = phase_space.positions
p_x, p_y, p_z = phase_space.momenta
phase_space.get_position(x) # Result: p_x
phase_space.get_momentum(p_y) # Result: y
```

Finally, the `PhaseSpace` object can be used to compute Hamilton's Equations

```python
hamiltons_equations = phase_space.hamiltons_equations()
print(hamiltons_equations)
```

Alternatively, we can get Hamilton's Equation's one coordinate or multiple coordinates at a time:

```python
x_equations = phase_space.hamiltons_equation(x, p_x)
```

We can use `sympy.dsolve` to solve the system of differential equations:

```python
from sympy import dsolve
# Define Initial Conditions
initial_conditions = {
    x.subs(t,0): Symbol("x_0"), # x(0) = x_0
    y.subs(t,0): Symbol("y_0"), # y(0) = y_0
    z.subs(t,0): Symbol("z_0"), # z(0) = z_0
    x.diff(t).subs(t,0): Symbol("v_x"), # x'(0) = v_x
    y.diff(t).subs(t,0): Symbol("v_y"), # y'(0) = v_y
    z.diff(t).subs(t,0): Symbol("v_{z0}"), # z'(0) = v_{z0}
}
# Use dsolve to solve the system
soln = dsolve(hamiltons_equations, ics=initial_conditions)
# Print the solution
print(soln)
```

## Example 2: A Lagrangian that depends on Accelerations

The `hamiltonian` package will compute Hamiltonians for Lagrangians that depend on higher order derivatives of the coordinates ($\ddot{x}, x^{(3)},...$). The `hamiltonian` package uses Ostrogradsky's method to compute these Hamiltonians.

Consider the following Lagrangian:

$$\mathcal{L} = -\frac{\epsilon m}{2\omega^2} \ddot{q}^2 + \frac{m}{2}\dot{q}^2-\frac{m\omega^2}{2}q^2$$

The Lagrangian above represents a simple harmonic oscillator which has been perturbed by a small dimensionless parameter $\epsilon$.

We can compute the Hamiltonian as in Example 1:

```python
epsilon, omega = symbol(["epsilon","omega"])
# Define Coordinate System
coordinate_system = hamiltonian.CoordinateSystem("t", ["q"])
t = coordinate_system.independent_variable
q = coordinate_system.get_coordinate("q")

# Define Lagrangian
L = - epsilon*m/(2*omega**2) * (q.diff(t,2))**2 + m/2 * (q.diff(t))**2 - m*omega**2/2 * q**2

# Compute Hamiltonian
phase_space = hamiltonian.hamiltonian(L, coordinate_system)
print(Eq("\mathcal{H}", phase_space.hamiltonian))

# Compute Hamilton's Equations
hamiltons_equations = phase_space.hamiltons_equations()
print(hamiltons_equations)
```

Note that because the Lagrangian depends on $\ddot{q}$, there are two "positions" in the resulting phase space ($q, q_1=\dot{q}$):

```python
print(phase_space.positions)
print(phase_space.momenta)
```
