# Elastic Leidenfrost

<img src="https://github.com/NCJCoombs/Elastic_Leidenfrost/blob/main/bulk_elastic_waves.png" width="800">

## Problem description

This code simulates the Leidenfrost effect for hydrogel spheres using the finite element method.

The user must first have an installation of the finite element library [oomph-lib](https://oomph-lib.github.io/oomph-lib/doc/html/).

The dimensional physical parameters for this problem are:

Parameter  | Symbol
------------- | -------------
Hydrogel Young's modulus | $E$ 
Hydrogel Poisson ratio | $\nu$ 
Hydrogel density | $\rho_s$ 
Vapour density | $\rho_v$ 
Vapour viscosity | $\mu$ 
Latent heat of evaporation | $L$
Substrate-hydrogel temperature difference | $\Delta T$ 
Vapour thermal conductivity | $k_v$

The dimensionless qunatities to be specified in the code, which may be passed as command line arguments, are:

Dimensionless parameter | Symbol | Formula | Command line argument
------------- | ------------- | ------------- | -------------
Leidenfrost source term | $\mathcal{E}$ |  $\left(\frac{\rho_s}{E}\right)^\frac{1}{2}\frac{k_v\Delta T}{\rho_v L a}$ | ```--conduction_source```
Pre-factor multiplying external traction | $\alpha$ | $\mu/[a(\rho_s E)^\frac{1}{2}]$ | ```--alpha```
Dimensionless gravitional forcing | $\bar{g}$ | $\rho_s a g/E$ | ```--g_non_dim```

Additional parameters, which have default values if not otherwise specified, are:

Parameter  | Command line argument
------------- | -------------
Initial height of the drop centre | ```--initial_height```
Initial drop speed | ```--initial_speed```
Boolean flag for use of solid pressure formaulation | ```--use_solid_pressure```
Boolean flag for use of incompressible formulation | ```--incompressible```

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
