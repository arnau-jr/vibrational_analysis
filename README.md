# Vibrational analysis package
## Description
This Fortran 90 code contains tools for the analysis and study of the vibrational dynamics of a molecule. Given an equilibrium configuration and force field this package allows:

- Evaluation of the potential energy of the molecule.
- Computation of the eigenvectors and eigenvalues of the normal mode representation of the molecule.
- Excitations of said normal modes to generate new configurations.
- Analyze the instantaneous configuration of a molecule from a MD simulation to obtain the Eckart frame, the kinetic energy decomposition and the normal mode energies of the molecule at that instant.

## Structure and compilation
The code is structured in three Fortran modules:
- Molecule: contains all the information of the inputted molecule (coordinates, velocities, masses, etc.)
- Force Field: contains all the information of the force field like number of bonds, angles, dihedrals and impropers; their types and coefficients; and all the necessary tools to evaluate the potential energy of the molecule (distance and angle computations, etc.). This module also contains the subroutines that compute the gradient and Hessian of the potential energy.
- Vibration: contains various subroutines for vibrational analysis. Most notably the computation of normal modes, Eckart frame, kinetic energy decomposition and excitation subroutines.

Vibration depends both on Molecule and Force Field but Molecule and Force Field are independent.

## Example usage
In the `example` directory there is a small script which takes a molecule with a force field as input, computes the normal modes and excites each of them to show the deformation associated to each of them. Included in the directory there is a nitromethane molecule files with a force field to use as inputs.

To compile you should place yourself in the `vibrational_analysis` and type `make all`. This will compile the modules. Then move to the `example` directory and type `make` (it may be necessary to perform this last make twice if the first one throws an error).

Executing `./test.x nitro_eq.xyz nitro_ff.dat` will then output each trajectory into the `results directory`.

## Force field structure and modification
######TODO

## Author and contact
Arnau Jurado Romero: arnau.jurado.romero@gmail.com
