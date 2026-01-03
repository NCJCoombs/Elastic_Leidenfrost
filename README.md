# Elastic Leidenfrost

<img src="https://github.com/NCJCoombs/Elastic_Leidenfrost/blob/main/bulk_elastic_waves.png" width="800">

## Problem description

This code simulates the Leidenfrost effect for hydrogel spheres using the finite element method.

The user must first have an installation of the finite element library [oomph-lib](https://oomph-lib.github.io/oomph-lib/doc/html/).

The dimensional physical parameters for this problem are:

Parameter  | Symbol | Default value
------------- | ------------- | -
Hydrogel Young's modulus | $E$ | $\pu{50 kPa}$
Hydrogel Poisson ratio | $\nu$ |$0.45$
Hydrogel density | $\rho_s$ |$\pu{1000 kg m^{-3}}$
Vapour density | $\rho_v$ | $\pu{0.5 kg m^{-3}}$
Vapour viscosity | $\mu$ | $\pu{2.0e-5 Pa s}$
Latent heat of evaporation | $L$ | $\pu{2.6e6 J kg^{-1}}$
Substrate-hydrogel temperature difference | $\Delta T$ | $\pu{115 K}$
Vapour thermal conductivity | $k_v$ | $\pu{0.03 W m^{-1} K^{-1}}$

The dimensionless qunatities to be specified in the code, which may be passed as command line arguments, are:

Dimensionless parameter | Symbol | Formula | Command line argument
------------- | ------------- | ------------- | -------------
Leidenfrost source term | $\mathcal{E}$ |  $\left(\frac{\rho_s}{E}\right)^\frac{1}{2}\frac{k_v\Delta T}{\rho_v L a}$ | ```--evap_number```
Pre-factor multiplying external traction | $\alpha$ | $\mu/[a(\rho_s E)^\frac{1}{2}]$ | ```--alpha```
Dimensionless gravitional forcing | $\bar{g}$ | $\rho_s a g/E$ | ```--gravity```

If not specified by the user, these parameters are calculated using the dimensional values. Additional parameters, which have default values if not otherwise specified, are:

Parameter  | Command line argument | Default value in code
------------- | ------------- | -
Initial height of the drop centre | ```--initial_height``` | $1.1$
Initial drop speed | ```--initial_speed``` | $\pu{1.4e-4}$
Boolean flag for use of solid pressure formaulation | ```--use_solid_pressure``` | $0$
Boolean flag for use of incompressible formulation | ```--incompressible``` | $0$

Other parameters, which control the tolerances for adaptive refinement and adaptive timestepping, for instance, cannot be altered using a command line flag, but can be changed within the driver itself.

## File description

### ```src/meshes```
```half_circle_sector_domain.h```, ```half_circle_sector_mesh.template.h```, ```half_circle_sector_mesh.template.cc```: Collectively defines the geometry of the semi-circular mesh in the undeformed configuration.
### ```user_src/axisym_cylindrical_solid```
```(refineable_)axisym_cylindrical_solid_elements.h```, ```(refineable_)axisym_cylindrical_solid_elements.cc```: Calculates the bulk residual and Jacobian contributions at the element level for non-linear elasticity in cylindrical coordinates with axial symmetry. The formulation used is the principle of virtual displacements; see [here](https://oomph-lib.github.io/oomph-lib/doc/solid/solid_theory/html/index.html) for more information. The refineable versions of these elements account for the presence of hanging nodes introduced by quadtree refinement.

```neo_hookean_constitutive.h```: Implements the strain energy density function, as well as its derivatives with respect to the strain invariants, for a Neo-Hookean material.

```solid_lubrication_elements.h```: Calculates the surface residual and Jacobian contributions at the element level that arise from the lubricating vapour film.
### ```user_drivers/elastic_leidenfrost```
```dynamic_EL_gnl_structured.cc```: The driver code for simulating a full isothermal or Leidenfrost impact.

```dynamic_EL_gnl_structured_steady.cc```: The driver code for obtaining stationary floating configurations by first simulating a low velocity impact and applying mass damping.

## Installation

The code has been written for and successfully compiled with [version 2.0.0](https://github.com/oomph-lib/oomph-lib/releases/tag/v2.0.0) of oomph-lib. The directories for the elements and driver codes mimic the directory structure of oomph-lib. Once these files have been placed in the relevant locations, and the name ```half_circle_sector_mesh``` has been appended to the file ```src/meshes/mesh_names.list```, the user may run ```./autogen.sh``` to generate the makefile for both drivers.

It is recommended that the user also installs the [mumps](https://mumps-solver.org/index.php) linear solver, rather than relying on oomph-lib's default, SuperLU.

## Output and example run

At the end of each computatonal timestep, the code will either produce a reduced output of the nodal position and vapour pressure at the hydrogel-air interface, or a full output which also contains the strain field within the bulk of the hydrogel. The frequency of a full output can be altered by changing the value of ```steps_per_bulk_doc``` in the code. The file ```time.dat``` also contains the following data which aid post-processing:

Value  | Symbol
------------- | -------------
Time | $t$
Gravitational potential energy of hydrogel | $E_g$
Kinetic energy of hydrogel | $E_k$
Elastic potential energy of hydrogel | $E_e$
Power due to surface forces | $\mathcal{P}_s$
Radial position of the 'neck' of the deformed interface | $l$
Average interface height from $r=0$ to $r=l$ | $\langle h\rangle$

As an example run, consider the following:

```./dynamic_EL_gnl_structured --initial_height 1.2 --initial_speed 0 --alpha 3e-7 --poisson 0.495 --gravity 0.0014715 --evap_number 6e-5 --use_solid_pressure 1 --incompressible 0```

This command corresponds to an impact simulation in which the centre of the hydrogel dropped from an initial height of 1.2 times the sphere radius. The solid pressure formulation for near-incompressible solids is used. However, true incompressibility (which also requires a Poisson ratio of 0.5) is not enforced.
