// Non-inline functions for axisymmetric solid mechanics elements

#include "axisym_cylindrical_solid_elements.h"
#include "axisym_cylindrical_solid_with_pressure_elements.h"

namespace oomph
{
  /// Static default value for the timescale ratio
  double AxisymmetricCylindricalPVDEquations::Default_lambda_sq_value = 1.0;

  /// Static default value for the damping parameters
  double AxisymmetricCylindricalPVDEquations::Default_eta_value = 0.0;

  /// Static default value for the timescale ratio
  double
    AxisymmetricCylindricalPVDWithPressureEquations::Default_lambda_sq_value =
      1.0;

  /// Static default value for the damping parameters
  double AxisymmetricCylindricalPVDWithPressureEquations::Default_eta_value =
    0.0;

  /// Return the residuals for the equations of solid mechanics
  void AxisymmetricCylindricalPVDEquations::
    fill_in_contribution_to_residuals_axisym_pvd(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian,
                                                 const unsigned& flag)
  {
    // Set the number of Lagrangian coordinates
    unsigned n_lagrangian = 2;
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Integers to store local equation number and local unknown
    int local_eqn = 0;
    int local_unknown = 0;

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidxi(n_node, n_lagrangian);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Get Rayleigh damping parameters
    const double eta_M = eta_mass();
    const double eta_K = eta_stiffness();

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    //----------------Compute the mass and stiffness matrices first--------//

    // Mass matrix, initialize all entries to zero. Size is 2*n_node as there
    // are two displacement components we care about
    DenseMatrix<double> M(2 * n_node, 2 * n_node, 0.0);
    // Stiffness matrix, initialize all entries to zero
    DenseMatrix<double> K(2 * n_node, 2 * n_node, 0.0);
    // Body force vector, initialize all entries to zero
    Vector<double> F(2 * n_node, 0.0);
    // Initialize derivatives of the stiffness matrix
    RankThreeTensor<double> dK_dX(2 * n_node, 2 * n_node, 2 * n_node, 0.0);

    // Get the timestepper from the first node
    TimeStepper* time_step_pt = node_pt(0)->time_stepper_pt();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);
      // Call the derivatives of the shape functions
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate the local Lagrangian coordinates, position components
      // and the derivatives of global position components
      // wrt lagrangian coordinates, as well as acceleration
      Vector<double> interpolated_xi(2, 0.0);
      Vector<double> interpolated_X(2, 0.0);
      DenseMatrix<double> interpolated_dXdxi(2, 2, 0.0);

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < 2; i++)
        {
          // Set the value of the lagrangian coordinate
          interpolated_xi[i] += lagrangian_position(l, i) * psi(l);
          // Set the value of the position component
          interpolated_X[i] += nodal_position(l, i) * psi(l);
          // Loop over Lagrangian derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            // Calculate dX[i]/dxi_{j}
            interpolated_dXdxi(i, j) += nodal_position(l, i) * dpsidxi(l, j);
          }
        }
      }

      // We are now in a position to calculate the undeformed metric tensor
      DenseMatrix<double> g(3);
      // r row
      g(0, 0) = 1.0;
      g(0, 1) = 0.0;
      g(0, 2) = 0.0;
      // z row
      g(1, 0) = 0.0;
      g(1, 1) = 1.0;
      g(1, 2) = 0.0;
      // phi row
      g(2, 0) = 0.0;
      g(2, 1) = 0.0;
      g(2, 2) = interpolated_xi[0] * interpolated_xi[0];

      // Now multiply the weight by the square-root of the undeformed metric
      // tensor r
      W *= sqrt(g(0, 0) * g(1, 1) * g(2, 2));

      // Now calculate the deformed metric tensor
      DenseMatrix<double> G(3);
      // r row
      G(0, 0) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 0) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 0);
      G(0, 1) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 1);
      G(0, 2) = 0.0;
      // z row
      G(1, 0) = G(0, 1);
      G(1, 1) = interpolated_dXdxi(0, 1) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 1) * interpolated_dXdxi(1, 1);
      G(1, 2) = 0.0;
      // phi row
      G(2, 0) = 0.0;
      G(2, 1) = 0.0;
      G(2, 2) = interpolated_X[0] * interpolated_X[0];

      // Now calculate the stress tensor from the constitutive law
      DenseMatrix<double> sigma(3, 3, 0.0);
      get_stress(g, G, sigma);

      // If we're calculating the Jacobian, will need derivative of stress
      // tensor w.r.t. the deformed metric tensor
      RankFourTensor<double> d_stress_dG(3, 3, 3, 3, 0.0);
      RankThreeTensor<double> dG_dR(n_node, 3, 3, 0.0);
      RankThreeTensor<double> dG_dZ(n_node, 3, 3, 0.0);

      if (flag == 1)
      {
        // Get the derivatives
        this->get_d_stress_dG_upper(g, G, sigma, d_stress_dG);

        // Loop over nodes
        for (unsigned m = 0; m < n_node; m++)
        {
          // Populate entries of dGdX matrices
          // R components
          dG_dR(m, 0, 0) = 2.0 * interpolated_dXdxi(0, 0) * dpsidxi(m, 0);
          dG_dR(m, 0, 1) = interpolated_dXdxi(0, 0) * dpsidxi(m, 1) +
                           interpolated_dXdxi(0, 1) * dpsidxi(m, 0);
          dG_dR(m, 1, 0) = dG_dR(m, 0, 1);
          dG_dR(m, 1, 1) = 2.0 * interpolated_dXdxi(0, 1) * dpsidxi(m, 1);
          dG_dR(m, 2, 2) = 2.0 * interpolated_X[0] * psi(m);

          // Z components
          dG_dZ(m, 0, 0) = 2.0 * interpolated_dXdxi(1, 0) * dpsidxi(m, 0);
          dG_dZ(m, 0, 1) = interpolated_dXdxi(1, 0) * dpsidxi(m, 1) +
                           interpolated_dXdxi(1, 1) * dpsidxi(m, 0);
          dG_dZ(m, 1, 0) = dG_dZ(m, 0, 1);
          dG_dZ(m, 1, 1) = 2.0 * interpolated_dXdxi(1, 1) * dpsidxi(m, 1);
        }
      }

      // Get body force at current time
      Vector<double> b(2, 0.0);
      this->body_force(interpolated_xi, b);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Forcing contributions
        F[l] -= psi(l) * b[0] * W;
        F[l + n_node] -= psi(l) * b[1] * W;

        // Second loop over test functions
        for (unsigned m = 0; m < n_node; m++)
        {
          //----RR COMPONENT---------------------//

          // Mass matrix contributions
          M(l, m) += psi(l) * lambda_sq * psi(m) * W;

          // Stiffness matrix contributions
          for (unsigned ii = 0; ii < 2; ii++)
          {
            for (unsigned jj = 0; jj < 2; jj++)
            {
              K(l, m) += sigma(ii, jj) * dpsidxi(l, jj) * dpsidxi(m, ii) * W;
            }
          }
          K(l, m) += sigma(2, 2) * psi(l) * psi(m) * W;

          //----RZ COMPONENT---------------------//

          // (None)

          //----ZR COMPONENT---------------------//

          // (None)

          //----ZZ COMPONENT---------------------//

          // Mass matrix contributions
          M(l + n_node, m + n_node) += psi(l) * lambda_sq * psi(m) * W;

          // Stiffness matrix contributions
          for (unsigned ii = 0; ii < 2; ii++)
          {
            for (unsigned jj = 0; jj < 2; jj++)
            {
              K(l + n_node, m + n_node) +=
                sigma(ii, jj) * dpsidxi(l, jj) * dpsidxi(m, ii) * W;
            }
          }

          // If we're doing Jacobian we'll need derivatives of the stiffness
          // matrix too
          if (flag)
          {
            // Derivative of stiffness matrix contributions
            for (unsigned k = 0; k < n_node; k++)
            {
              // Loop over deformed metric tensor components
              for (unsigned alpha = 0; alpha < 3; alpha++)
              {
                for (unsigned beta = alpha; beta < 3; beta++)
                {
                  double pre_factor = 1.0;
                  if (alpha != beta) pre_factor *= 2.0;

                  double term =
                    pre_factor * (d_stress_dG(0, 0, alpha, beta) *
                                    dpsidxi(l, 0) * dpsidxi(m, 0) +
                                  d_stress_dG(0, 1, alpha, beta) *
                                    (dpsidxi(l, 0) * dpsidxi(m, 1) +
                                     dpsidxi(l, 1) * dpsidxi(m, 0)) +
                                  d_stress_dG(1, 1, alpha, beta) *
                                    dpsidxi(l, 1) * dpsidxi(m, 1));

                  dK_dX(k, l, m) += term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l, m) += term * dG_dZ(k, alpha, beta) * W;

                  dK_dX(k, l + n_node, m + n_node) +=
                    term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l + n_node, m + n_node) +=
                    term * dG_dZ(k, alpha, beta) * W;

                  // Due to axisymmetry
                  term = pre_factor * d_stress_dG(2, 2, alpha, beta) * psi(l) *
                         psi(m);

                  dK_dX(k, l, m) += term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l, m) += term * dG_dZ(k, alpha, beta) * W;
                }
              }
            }
          }
        }
      }
    } // End of loop over integration points

    // Rayleigh damping matrix, composed of the mass and stiffness matrices
    DenseMatrix<double> C(2 * n_node, 2 * n_node, 0.0);
    for (unsigned l = 0; l < 2 * n_node; l++)
    {
      for (unsigned m = 0; m < 2 * n_node; m++)
      {
        C(l, m) = eta_M * M(l, m) + eta_K * K(l, m);
      }
    }

    // Re-do node loops and assemble residuals/Jacobian
    // First node loop
    for (unsigned l = 0; l < n_node; l++)
    {
      //---------------R residual------------------------------//

      // The middle 0 here denotes the zeroth position type (assuming
      // Lagrange-type elements)
      local_eqn = position_local_eqn(l, 0, 0);

      // If it's not a boundary condition
      if (local_eqn >= 0)
      {
        // Do force contribution
        residuals[local_eqn] += F[l];

        // Second node loop
        for (unsigned m = 0; m < n_node; m++)
        {
          // Add the residual contribution RR
          residuals[local_eqn] += M(l, m) * dnodal_position_dt(m, 2, 0) +
                                  C(l, m) * dnodal_position_dt(m, 1, 0) +
                                  K(l, m) * nodal_position(m, 0);

          // If we're computing the jacobian too
          if (flag)
          {
            //--------RR contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 0);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                M(l, m) * time_step_pt->weight(2, 0) +
                C(l, m) * time_step_pt->weight(1, 0) + K(l, m);

              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m, l, k) *
                  (eta_K * dnodal_position_dt(k, 0) + nodal_position(k, 0));
              }
            }

            //--------RZ contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 1);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m + n_node, l, k) *
                  (eta_K * dnodal_position_dt(k, 0) + nodal_position(k, 0));
              }
            }
          }
        }
      }

      //---------------Z residual------------------------------//

      // The middle 0 here denotes the zeroth position type (assuming
      // Lagrange-type elements)
      local_eqn = position_local_eqn(l, 0, 1);

      // If it's not a boundary condition
      if (local_eqn >= 0)
      {
        // Do force contribution
        residuals[local_eqn] += F[l + n_node];

        // Second node loop
        for (unsigned m = 0; m < n_node; m++)
        {
          // Add the residual contribution ZZ
          residuals[local_eqn] +=
            M(l + n_node, m + n_node) * dnodal_position_dt(m, 2, 1) +
            C(l + n_node, m + n_node) * dnodal_position_dt(m, 1, 1) +
            K(l + n_node, m + n_node) * nodal_position(m, 1);

          // If we're computing the jacobian too
          if (flag)
          {
            //--------ZR contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 0);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m, l + n_node, k + n_node) *
                  (eta_K * dnodal_position_dt(k, 1) + nodal_position(k, 1));
              }
            }

            //--------ZZ contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 1);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                M(l + n_node, m + n_node) * time_step_pt->weight(2, 0) +
                C(l + n_node, m + n_node) * time_step_pt->weight(1, 0) +
                K(l + n_node, m + n_node);

              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m + n_node, l + n_node, k + n_node) *
                  (eta_K * dnodal_position_dt(k, 1) + nodal_position(k, 1));
              }
            }
          }
        }
      }
    }
  }

  void AxisymmetricCylindricalPVDWithPressureEquations::
    fill_in_contribution_to_residuals_axisym_pvd_with_pressure(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // Set the number of Lagrangian coordinates
    unsigned n_lagrangian = 2;
    // Find out how many nodes there are
    unsigned n_node = nnode();
    // Find out how many pressure dofs there are
    unsigned n_solid_pres = nsolid_pres();

    // Integers to store local equation number and local unknown
    int local_eqn = 0;
    int local_unknown = 0;

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsidxi(n_node, n_lagrangian);

    // Set up memory for the pressure shape functions
    Shape psisp(n_solid_pres);

    // Timescale ratio (non-dim density)
    double lambda_sq = this->lambda_sq();

    // Get Rayleigh damping parameters
    const double eta_M = eta_mass();
    const double eta_K = eta_stiffness();

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    //----------------Compute the mass and stiffness matrices first--------//

    // Mass matrix, initialize all entries to zero. Size is 2*n_node as there
    // are two displacement components we care about
    DenseMatrix<double> M(2 * n_node, 2 * n_node, 0.0);
    // Stiffness matrix, initialize all entries to zero
    DenseMatrix<double> K(2 * n_node, 2 * n_node, 0.0);
    // Body force vector, initialize all entries to zero
    Vector<double> F(2 * n_node, 0.0);
    // Initialize derivatives of the stiffness matrix
    RankThreeTensor<double> dK_dX(2 * n_node, 2 * n_node, 2 * n_node, 0.0);

    // Get the timestepper from the first node
    TimeStepper* time_step_pt = node_pt(0)->time_stepper_pt();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);
      // Call the derivatives of the shape functions
      double J = dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
      // Call the pressure shape functions
      solid_pshape_at_knot(ipt, psisp);
      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Calculate the local Lagrangian coordinates, position components
      // and the derivatives of global position components
      // wrt lagrangian coordinates, as well as acceleration
      Vector<double> interpolated_xi(2, 0.0);
      Vector<double> interpolated_X(2, 0.0);
      DenseMatrix<double> interpolated_dXdxi(2, 2, 0.0);
      double interpolated_solid_p = 0.0;

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over displacement components (deformed position)
        for (unsigned i = 0; i < 2; i++)
        {
          // Set the value of the lagrangian coordinate
          interpolated_xi[i] += lagrangian_position(l, i) * psi(l);
          // Set the value of the position component
          interpolated_X[i] += nodal_position(l, i) * psi(l);
          // Loop over Lagrangian derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            // Calculate dX[i]/dxi_{j}
            interpolated_dXdxi(i, j) += nodal_position(l, i) * dpsidxi(l, j);
          }
        }
      }

      // Calculate the local internal pressure
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        interpolated_solid_p += solid_p(l) * psisp(l);
      }

      // We are now in a position to calculate the undeformed metric tensor
      DenseMatrix<double> g(3);
      // r row
      g(0, 0) = 1.0;
      g(0, 1) = 0.0;
      g(0, 2) = 0.0;
      // z row
      g(1, 0) = 0.0;
      g(1, 1) = 1.0;
      g(1, 2) = 0.0;
      // phi row
      g(2, 0) = 0.0;
      g(2, 1) = 0.0;
      g(2, 2) = interpolated_xi[0] * interpolated_xi[0];

      // Now multiply the weight by the square-root of the undeformed metric
      // tensor r
      double detg = g(0, 0) * g(1, 1) * g(2, 2);
      W *= sqrt(detg);

      // Now calculate the deformed metric tensor
      DenseMatrix<double> G(3);
      // r row
      G(0, 0) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 0) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 0);
      G(0, 1) = interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 1);
      G(0, 2) = 0.0;
      // z row
      G(1, 0) = G(0, 1);
      G(1, 1) = interpolated_dXdxi(0, 1) * interpolated_dXdxi(0, 1) +
                interpolated_dXdxi(1, 1) * interpolated_dXdxi(1, 1);
      G(1, 2) = 0.0;
      // phi row
      G(2, 0) = 0.0;
      G(2, 1) = 0.0;
      G(2, 2) = interpolated_X[0] * interpolated_X[0];

      // Now calculate the deviatoric stress tensor from the constitutive law
      DenseMatrix<double> sigma_dev(3), Gup(3);
      double detG = 0.0, gen_dil = 0.0, inv_kappa = 0.0;
      // If it's incompressible call one form of the constitutive law
      if (Incompressible)
      {
        get_stress(g, G, sigma_dev, Gup, detG);
      }
      // Otherwise call another form
      else
      {
        get_stress(g, G, sigma_dev, Gup, gen_dil, inv_kappa);
      }

      // Build the stress tensor up from its pressure and deviatoric
      // components
      DenseMatrix<double> sigma(3, 3, 0.0);
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          sigma(i, j) =
            -1.0 * interpolated_solid_p * Gup(i, j) + sigma_dev(i, j);
        }
      }

      // If we're calculating the Jacobian, will need derivative of stress
      // tensor w.r.t. the deformed metric tensor
      RankFourTensor<double> d_stress_dG(3, 3, 3, 3, 0.0);
      RankThreeTensor<double> dG_dR(n_node, 3, 3, 0.0);
      RankThreeTensor<double> dG_dZ(n_node, 3, 3, 0.0);
      DenseMatrix<double> d_detG_dG(3, 3, 0.0);
      DenseMatrix<double> d_gen_dil_dG(3, 3, 0.0);

      if (flag == 1)
      {
        // If incompressible, call the incompressible form
        if (Incompressible)
        {
          this->get_d_stress_dG_upper(
            g, G, sigma, detG, interpolated_solid_p, d_stress_dG, d_detG_dG);
        }
        else
        // Otherwise call the near-incompressible form
        {
          this->get_d_stress_dG_upper(g,
                                      G,
                                      sigma,
                                      gen_dil,
                                      inv_kappa,
                                      interpolated_solid_p,
                                      d_stress_dG,
                                      d_gen_dil_dG);
        }

        // Loop over nodes
        for (unsigned m = 0; m < n_node; m++)
        {
          // Populate entries of dGdX matrices
          // R components
          dG_dR(m, 0, 0) = 2.0 * interpolated_dXdxi(0, 0) * dpsidxi(m, 0);
          dG_dR(m, 0, 1) = interpolated_dXdxi(0, 0) * dpsidxi(m, 1) +
                           interpolated_dXdxi(0, 1) * dpsidxi(m, 0);
          dG_dR(m, 1, 0) = dG_dR(m, 0, 1);
          dG_dR(m, 1, 1) = 2.0 * interpolated_dXdxi(0, 1) * dpsidxi(m, 1);
          dG_dR(m, 2, 2) = 2.0 * interpolated_X[0] * psi(m);

          // Z components
          dG_dZ(m, 0, 0) = 2.0 * interpolated_dXdxi(1, 0) * dpsidxi(m, 0);
          dG_dZ(m, 0, 1) = interpolated_dXdxi(1, 0) * dpsidxi(m, 1) +
                           interpolated_dXdxi(1, 1) * dpsidxi(m, 0);
          dG_dZ(m, 1, 0) = dG_dZ(m, 0, 1);
          dG_dZ(m, 1, 1) = 2.0 * interpolated_dXdxi(1, 1) * dpsidxi(m, 1);
        }
      }

      // Get body force at current time
      Vector<double> b(2, 0.0);
      this->body_force(interpolated_xi, b);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Forcing contributions
        F[l] -= psi(l) * b[0] * W;
        F[l + n_node] -= psi(l) * b[1] * W;

        // Second loop over test functions
        for (unsigned m = 0; m < n_node; m++)
        {
          //----RR COMPONENT---------------------//

          // Mass matrix contributions
          M(l, m) += psi(l) * lambda_sq * psi(m) * W;

          // Stiffness matrix contributions
          for (unsigned ii = 0; ii < 2; ii++)
          {
            for (unsigned jj = 0; jj < 2; jj++)
            {
              K(l, m) += sigma(ii, jj) * dpsidxi(l, jj) * dpsidxi(m, ii) * W;
            }
          }
          K(l, m) += sigma(2, 2) * psi(l) * psi(m) * W;

          //----RZ COMPONENT---------------------//

          // (None)

          //----ZR COMPONENT---------------------//

          // (None)

          //----ZZ COMPONENT---------------------//

          // Mass matrix contributions
          M(l + n_node, m + n_node) += psi(l) * lambda_sq * psi(m) * W;

          // Stiffness matrix contributions
          for (unsigned ii = 0; ii < 2; ii++)
          {
            for (unsigned jj = 0; jj < 2; jj++)
            {
              K(l + n_node, m + n_node) +=
                sigma(ii, jj) * dpsidxi(l, jj) * dpsidxi(m, ii) * W;
            }
          }

          // If we're doing Jacobian we'll need derivatives of the stiffness
          // matrix too
          if (flag)
          {
            // Derivative of stiffness matrix contributions
            for (unsigned k = 0; k < n_node; k++)
            {
              // Loop over deformed metric tensor components
              for (unsigned alpha = 0; alpha < 3; alpha++)
              {
                for (unsigned beta = alpha; beta < 3; beta++)
                {
                  double pre_factor = 1.0;
                  if (alpha != beta) pre_factor *= 2.0;

                  double term =
                    pre_factor * (d_stress_dG(0, 0, alpha, beta) *
                                    dpsidxi(l, 0) * dpsidxi(m, 0) +
                                  d_stress_dG(0, 1, alpha, beta) *
                                    (dpsidxi(l, 0) * dpsidxi(m, 1) +
                                     dpsidxi(l, 1) * dpsidxi(m, 0)) +
                                  d_stress_dG(1, 1, alpha, beta) *
                                    dpsidxi(l, 1) * dpsidxi(m, 1));

                  dK_dX(k, l, m) += term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l, m) += term * dG_dZ(k, alpha, beta) * W;

                  dK_dX(k, l + n_node, m + n_node) +=
                    term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l + n_node, m + n_node) +=
                    term * dG_dZ(k, alpha, beta) * W;

                  // Due to axisymmetry
                  term = pre_factor * d_stress_dG(2, 2, alpha, beta) * psi(l) *
                         psi(m);

                  dK_dX(k, l, m) += term * dG_dR(k, alpha, beta) * W;
                  dK_dX(k + n_node, l, m) += term * dG_dZ(k, alpha, beta) * W;
                }
              }
            }

            // Loop over directions
            for (unsigned i = 0; i < 2; i++)
            {
              // The local eqn for displacement
              local_eqn = position_local_eqn(l, 0, i);

              // If not a boundary condition
              if (local_eqn >= 0)
              {
                // Loop over solid pressure degrees of freedom
                for (unsigned ll = 0; ll < n_solid_pres; ll++)
                {
                  // Get the local solid unknown
                  local_unknown = solid_p_local_eqn(ll);

                  // If not a boundary condition
                  if (local_unknown >= 0)
                  {
                    // Loop over directions
                    for (unsigned ii = 0; ii < 2; ii++)
                    {
                      for (unsigned jj = 0; jj < 2; jj++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          psisp(ll) * Gup(ii, jj) * dpsidxi(l, jj) *
                          dpsidxi(m, ii) * W * nodal_position(m, i);
                      }
                    }

                    // If it's the radial local eqn, then there is also the
                    // phi-phi term
                    if (i == 0)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        psisp(ll) * Gup(2, 2) * psi(l) * psi(m) * W *
                        nodal_position(m, i);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Loop over the pressure degrees of freedom
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        local_eqn = solid_p_local_eqn(l);
        // If it's not a bondary condition
        if (local_eqn >= 0)
        {
          // For true incompressibility we need the ratio of determinants of
          // the metric tensors to be exactly 1.0
          if (Incompressible)
          {
            residuals[local_eqn] += (detG / detg - 1.0) * psisp(l) * W;

            // Add in Jacobian terms
            if (flag)
            {
              // Loop over nodes
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over directions
                for (unsigned i = 0; i < 2; i++)
                {
                  local_unknown = position_local_eqn(ll, 0, i);

                  // If not a boundary condition
                  if (local_unknown >= 0)
                  {
                    // Loop over deformed metric tensor components
                    for (unsigned alpha = 0; alpha < 3; alpha++)
                    {
                      for (unsigned beta = alpha; beta < 3; beta++)
                      {
                        double pre_factor = 1.0;
                        if (alpha != beta) pre_factor *= 2.0;

                        if (i == 0)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            pre_factor * (1.0 / detg) * d_detG_dG(alpha, beta) *
                            dG_dR(ll, alpha, beta) * psisp(l) * W;
                        }
                        else
                        {
                          jacobian(local_eqn, local_unknown) +=
                            pre_factor * (1.0 / detg) * d_detG_dG(alpha, beta) *
                            dG_dZ(ll, alpha, beta) * psisp(l) * W;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          else
          {
            // Otherwise the pressure must be that calculated by the
            // constitutive law
            residuals[local_eqn] +=
              (inv_kappa * interpolated_solid_p + gen_dil) * psisp(l) * W;

            // Add in the jacobian terms
            if (flag)
            {
              // Loop over the pressure nodes again
              for (unsigned ll = 0; ll < n_solid_pres; ll++)
              {
                local_unknown = solid_p_local_eqn(ll);
                // If not a boundary condition
                if (local_unknown >= 0)
                {
                  jacobian(local_eqn, local_unknown) +=
                    inv_kappa * psisp(ll) * psisp(l) * W;
                }
              }

              // Loop over nodes
              for (unsigned ll = 0; ll < n_node; ll++)
              {
                // Loop over directions
                for (unsigned i = 0; i < 2; i++)
                {
                  local_unknown = position_local_eqn(ll, 0, i);

                  // If not a boundary condition
                  if (local_unknown >= 0)
                  {
                    // Loop over deformed metric tensor components
                    for (unsigned alpha = 0; alpha < 3; alpha++)
                    {
                      for (unsigned beta = alpha; beta < 3; beta++)
                      {
                        double pre_factor = 1.0;
                        if (alpha != beta) pre_factor *= 2.0;

                        if (i == 0)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            pre_factor * d_gen_dil_dG(alpha, beta) *
                            dG_dR(ll, alpha, beta) * psisp(l) * W;
                        }
                        else
                        {
                          jacobian(local_eqn, local_unknown) +=
                            pre_factor * d_gen_dil_dG(alpha, beta) *
                            dG_dZ(ll, alpha, beta) * psisp(l) * W;
                        }
                      }
                    }
                  }
                }
              }
            } // End of jacobian terms
          } // End of else
        } // End of if not boundary condition
      } // End of pressure dof loop
    } // End of loop over integration points

    // Rayleigh damping matrix, composed of the mass and stiffness matrices
    DenseMatrix<double> C(2 * n_node, 2 * n_node, 0.0);
    for (unsigned l = 0; l < 2 * n_node; l++)
    {
      for (unsigned m = 0; m < 2 * n_node; m++)
      {
        C(l, m) = eta_M * M(l, m) + eta_K * K(l, m);
      }
    }

    // Re-do node loops and assemble residuals/Jacobian
    // First node loop
    for (unsigned l = 0; l < n_node; l++)
    {
      //---------------R residual------------------------------//

      // The middle 0 here denotes the zeroth position type (assuming
      // Lagrange-type elements)
      local_eqn = position_local_eqn(l, 0, 0);

      // If it's not a boundary condition
      if (local_eqn >= 0)
      {
        // Do force contribution
        residuals[local_eqn] += F[l];

        // Second node loop
        for (unsigned m = 0; m < n_node; m++)
        {
          // Add the residual contribution RR
          residuals[local_eqn] += M(l, m) * dnodal_position_dt(m, 2, 0) +
                                  C(l, m) * dnodal_position_dt(m, 1, 0) +
                                  K(l, m) * nodal_position(m, 0);

          // If we're computing the jacobian too
          if (flag)
          {
            //--------RR contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 0);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                M(l, m) * time_step_pt->weight(2, 0) +
                C(l, m) * time_step_pt->weight(1, 0) + K(l, m);

              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m, l, k) *
                  (eta_K * dnodal_position_dt(k, 0) + nodal_position(k, 0));
              }
            }

            //--------RZ contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 1);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m + n_node, l, k) *
                  (eta_K * dnodal_position_dt(k, 0) + nodal_position(k, 0));
              }
            }
          }
        }
      }

      //---------------Z residual------------------------------//

      // The middle 0 here denotes the zeroth position type (assuming
      // Lagrange-type elements)
      local_eqn = position_local_eqn(l, 0, 1);

      // If it's not a boundary condition
      if (local_eqn >= 0)
      {
        // Do force contribution
        residuals[local_eqn] += F[l + n_node];

        // Second node loop
        for (unsigned m = 0; m < n_node; m++)
        {
          // Add the residual contribution ZZ
          residuals[local_eqn] +=
            M(l + n_node, m + n_node) * dnodal_position_dt(m, 2, 1) +
            C(l + n_node, m + n_node) * dnodal_position_dt(m, 1, 1) +
            K(l + n_node, m + n_node) * nodal_position(m, 1);

          // If we're computing the jacobian too
          if (flag)
          {
            //--------ZR contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 0);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m, l + n_node, k + n_node) *
                  (eta_K * dnodal_position_dt(k, 1) + nodal_position(k, 1));
              }
            }

            //--------ZZ contribution----------------------------//
            // Get the local unknown
            local_unknown = position_local_eqn(m, 0, 1);

            // if not a boundary condition
            if (local_unknown >= 0)
            {
              jacobian(local_eqn, local_unknown) +=
                M(l + n_node, m + n_node) * time_step_pt->weight(2, 0) +
                C(l + n_node, m + n_node) * time_step_pt->weight(1, 0) +
                K(l + n_node, m + n_node);

              for (unsigned k = 0; k < n_node; k++)
              {
                jacobian(local_eqn, local_unknown) +=
                  dK_dX(m + n_node, l + n_node, k + n_node) *
                  (eta_K * dnodal_position_dt(k, 1) + nodal_position(k, 1));
              }
            }
          }
        }
      }
    }
  }

  //====================================================================
  /// Data for the number of Variables at each node
  //====================================================================
  const unsigned
    QAxisymmetricCylindricalPVDWithPressureElement::Initial_Nvalue[9] = {
      1, 0, 1, 0, 0, 0, 1, 0, 1};

  //==========================================================================
  /// Conversion from pressure dof to Node number at which pressure is stored
  //==========================================================================
  const unsigned QAxisymmetricCylindricalPVDWithPressureElement::Pconv[4] = {
    0, 2, 6, 8};

  //=======================================================================
  /// Data for the number of variables at each node
  //=======================================================================
  const unsigned TAxisymCylindricalPVDWithPressureElement::Initial_Nvalue[6] = {
    1, 1, 1, 0, 0, 0};

  //=======================================================================
  /// Data for the pressure conversion array
  //=======================================================================
  const unsigned TAxisymCylindricalPVDWithPressureElement::Pconv[3] = {0, 1, 2};

  //==================================================================
  /// Solid pressure shape function evaluated at integration point
  //==================================================================
  void AxisymmetricCylindricalPVDWithPressureEquations::solid_pshape_at_knot(
    const unsigned& ipt, Shape& psi) const
  {
    // Storage for local coordinates of the integration point
    Vector<double> s(2);
    // Set the local coordinates
    for (unsigned i = 0; i < 2; i++)
    {
      s[i] = this->integral_pt()->knot(ipt, i);
    }
    // Get the shape function
    solid_pshape(s, psi);
  }

  // Build required templates
  template class TAxisymCylindricalPVDElement<3>;
  template class QAxisymmetricCylindricalPVDElement<3>;

} // namespace oomph
