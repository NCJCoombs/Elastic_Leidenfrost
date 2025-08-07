// Driver code for an elastic solid Leidenfrost impact problem with geometric
// non-linearity

// Generic routines
#include "generic.h"

// Axisymmetric non-linear solid equations in cylindrical coords
#include "axisym_cylindrical_solid.h"

// Constitutive laws
#include "constitutive.h"

// The mesh
#include "meshes/half_circle_sector_mesh.h"
#include "meshes/one_d_mesh.h"

using namespace std;
using namespace oomph;

//=======start_namespace==========================================
/// Global physical variables
//================================================================
namespace GlobalPhysicalVariables
{
  ///--------------Physical parameters (SI units)----------///

  // Sphere radius
  double R = 7.0e-3;

  // Elastic modulus
  double E = 5.0e4;

  // Poisson ratio
  double nu = 0.45; // 0.495;

  // Temperature difference (substrate-solid)
  double delta_T = 115;

  // Shear viscosity of vapour
  double mu_v = 2.0e-5;

  // Thermal conductivity of vapour
  double k_v = 3.0e-2;

  // Latent heat
  double L = 2.6e6;

  // Vapour density
  double rho_v = 5.0e-1;

  // Hydrogel density
  double rho_s = 1.0e3;

  // Gravitational accel.
  double g = 9.81;

  // Initial vertical velocity
  double W = 1.0e-3;

  ///------------------------------------------------------///

  /// Non-dimensional parameters
  // Conduction source term
  double conduction_source =
    (k_v * delta_T / (L * rho_v * R)) * pow(rho_s / E, 0.5);

  // Pressure scaling parameter
  double pressure_scaling_param = mu_v / (R * pow(rho_s * E, 0.5));

  // Scaling factor in linear elasticity eqns (1 due to scaling)
  double lambda_sq = 1.0;

  // Non-dim gravitational forcing
  double g_non_dim = rho_s * R * g / E;

  // Non-dim initial velocity
  double W_non_dim = W / pow(E / rho_s, 0.5);

  // Forcing function due to gravity
  void get_body_force(const double& time,
                      const Vector<double>& x,
                      Vector<double>& b)
  {
    b[0] = 0.0;
    b[1] = -g_non_dim;
  }

  // Artificial repulsion
  void get_artificial_repulsion(const double& time, const double& h, double& f)
  {
    double epsilon = 1.0e6;
    double alpha = 2.0;
    double h_thres = 0.0001;
    f = epsilon * pow(h_thres / h, alpha);
  }

  // Initial height of the centre of the solid
  double h_initial = 1.005;

} // namespace GlobalPhysicalVariables

//=======start_namespace==========================================
/// Global simulation settings
//================================================================
namespace GlobalSimSettings
{
  /// Newmark parameters (0.5 for non-dissipative)
  double Newmark_beta1 = 0.5025;
  double Newmark_beta2 = 0.505;

  /// Flag for using a pressure formulation
  unsigned use_pressure_formulation = 0;

  /// If using a pressure formulation, is the solid incompressible
  unsigned incompressible = 0;

  /// Result Folder
  string result_folder = "RESLT";

  /// Height at which the simulation starts (not necessarily the initial solid
  /// height, initial condition is set accordingly)
  double h_initial_sim = 1.1;

  // Height at which the adaptive timestepping is turned on
  double height_impact = 0.25;

  // Height at which the adaptive remeshing is turned on
  double height_impact_var = 0.05;

  // Maximum timestep
  double max_dt = 0.001;

  // Targets for spatial adaptivity
  // Max error (elements refined if error is above this)
  double max_permitted_error = 5.0e-5;
  // Min error (elements unrefined if error is below this)
  double min_permitted_error = 1.0e-5;

  // Max refinement level
  unsigned max_refinement_level = 13;

  // Min refinement level
  unsigned min_refinement_level = 6;

  // Lubrication equation is solved  where
  // the z component of the normal to the drop is less than -film_cutoff
  double film_cutoff = 0.1;

  // Bool for use adaptive timestepping (initially false)
  bool use_adaptive_timestepping = false;

  // Temporal error tolerance
  double epsilon_t = 1.0e-4;

  // Pressure variation tolerance (set in timestepping loop)
  double epsilon_p = 1.0e3;

  // Number of computational timesteps per interface doc
  unsigned steps_per_interface_doc = 1;

  // Number of computational timesteps per bulk doc
  unsigned steps_per_bulk_doc = 2;

  // Number of steps per adapt
  unsigned steps_per_adapt = 3;

  // Boolean for whether base mesh has been created
  bool base_mesh_created = false;

} // namespace GlobalSimSettings

//===VariableNewmark====================================================
/// Newmark scheme with parameters that can be varied
//======================================================================
template<unsigned NSTEPS>
class VariableNewmark : public Newmark<NSTEPS>
{
public:
  // Constructor
  VariableNewmark()
  {
    this->Type = "Variable Newmark";
    this->Beta1 = GlobalSimSettings::Newmark_beta1;
    this->Beta2 = GlobalSimSettings::Newmark_beta2;
  }

  /// Broken copy constructor
  VariableNewmark(const VariableNewmark&) = delete;

  /// Broken assignment operator
  void operator=(const VariableNewmark&) = delete;

  /// Access function for setting the Newmark parameters
  void set_newmark_parameters(const double& beta1, const double& beta2)
  {
    this->Beta1 = beta1;
    this->Beta2 = beta2;
  }
};

//======================start_mesh================================
/// Elastic quarter circle sector mesh: We "upgrade"
/// the RefineableHalfCircleSectorMesh to become an
/// SolidMesh and equate the Eulerian and Lagrangian coordinates,
/// thus making the domain represented by the mesh the stress-free
/// configuration.
//================================================================
template<class ELEMENT>
class ElasticRefineableHalfCircleSectorMesh
  : public virtual RefineableHalfCircleSectorMesh<ELEMENT>,
    public virtual SolidMesh
{
public:
  /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
  /// ones so that the initial configuration is the stress-free one.
  ElasticRefineableHalfCircleSectorMesh<ELEMENT>(
    GeomObject* wall_pt,
    double x_centr,
    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
    : RefineableHalfCircleSectorMesh<ELEMENT>(wall_pt, x_centr, time_stepper_pt)
  {
    // Make the current configuration the undeformed one by
    // setting the nodal Lagrangian coordinates to their current
    // Eulerian ones
    set_lagrangian_nodal_coordinates();

    // Build "undeformed" domain: This is a "deep" copy of the
    // Domain that we used to create set the Eulerian coordinates
    // in the initial mesh -- the original domain (accessible via
    // the private member data Domain_pt) will be used to update
    // the position of the boundary nodes; the copy that we're
    // creating here will be used to determine the Lagrangian coordinates
    // of any newly created SolidNodes during mesh refinement
    Undeformed_domain_pt = new HalfCircleSectorDomain(wall_pt, x_centr);

    // Loop over all elements and set the undeformed macro element pointer
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to full element type
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));

      // Set pointer to macro element so the curvlinear boundaries
      // of the undeformed mesh/domain get picked up during adaptive
      // mesh refinement
      el_pt->set_undeformed_macro_elem_pt(
        Undeformed_domain_pt->macro_element_pt(e));

      // Use MacroElement representation for
      // Lagrangian coordinates of newly created
      // nodes
      el_pt->enable_use_of_undeformed_macro_element_for_new_lagrangian_coords();
    }
  }

  /// Destructor: Kill "undeformed" Domain
  virtual ~ElasticRefineableHalfCircleSectorMesh<ELEMENT>()
  {
    delete Undeformed_domain_pt;
  }


private:
  /// Pointer to "undeformed" Domain -- used to determine the
  /// Lagrangian coordinates of any newly created SolidNodes during
  /// Mesh refinement
  Domain* Undeformed_domain_pt;
};


//==============start_problem=========================================
/// Unstructured solid problem
//====================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
class ElasticLeidenfrostProblem : public Problem
{
public:
  /// Constructor:
  ElasticLeidenfrostProblem();

  /// Destructor
  ~ElasticLeidenfrostProblem()
  {
    delete outer_boundary_circle_pt;
    delete Bulk_mesh_pt;
    delete error_estimator_pt;
    delete Constitutive_law_pt;
    delete this->time_stepper_pt();
  }

  /// Estimate of temporal error
  double custom_global_temporal_error_norm();

  /// Custom wrapper for adaptive timestepping
  void custom_adaptive_unsteady_newton_solve(const double& dt_desired,
                                             const double& epsilon,
                                             double& dt_new,
                                             bool& success);

  /// Actions before adapt: Wipe the mesh of boundary elements
  void actions_before_adapt()
  {
    // New elements should not use the macro-element description if the base
    // mesh has already been created
    if (GlobalSimSettings::base_mesh_created)
    {
      // Switch to using FE interpolation for the creation of new nodes, as
      // opposed to the undeformed macro element representation
      unsigned n_el = Bulk_mesh_pt->nelement();
      for (unsigned e = 0; e < n_el; e++)
      {
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
        el_pt
          ->disable_use_of_undeformed_macro_element_for_new_lagrangian_coords();
      }
    }

    // Kill the  elements and wipe surface mesh
    delete_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();
  }

  /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
  void actions_after_adapt()
  {
    // Create the boundary elements
    create_boundary_elements();

    // Rebuild the Problem's global mesh from its various sub-meshes
    this->rebuild_global_mesh();

    // Setup the problem again -- remember that fluid mesh has been
    // completely rebuilt and its element's don't have any
    // pointers to Re etc.
    complete_problem_setup();

    // Loop over elements and pin solid redundant solid pressures (if any)
    unsigned n_el = Bulk_mesh_pt->nelement();
    for (unsigned e = 0; e < n_el; e++)
    {
      // Cast to element with pressure dofs
      RefineableQAxisymmetricCylindricalPVDWithPressureElement* el_pt =
        dynamic_cast<RefineableQAxisymmetricCylindricalPVDWithPressureElement*>(
          Bulk_mesh_pt->element_pt(e));

      // Cast will only succeed if we're using the pressure formulation
      if (el_pt != 0)
      {
        el_pt->pin_elemental_redundant_nodal_solid_pressures();
      }
    }
  }

  /// Set the problem's initial condition
  void set_initial_condition();

  /// Doc the solution
  void doc_solution();

  // Do we want to enable refinement in the bulk mesh
  void toggle_bulk_mesh_refinement(double& h_min)
  {
    // If past the threshold, turn on adaptive timestepping
    if (h_min < GlobalSimSettings::height_impact)
    {
      GlobalSimSettings::use_adaptive_timestepping = true;
    }
    else
    {
      GlobalSimSettings::use_adaptive_timestepping = false;
    }

    // If past the threshold, turn on adaptive remeshing
    if (h_min < GlobalSimSettings::height_impact_var)
    {
      Bulk_mesh_pt->enable_adaptation();
    }
    else
    {
      Bulk_mesh_pt->disable_adaptation();
    }

    // Get the contact radius and average height
    double l, h_avg;
    compute_contact_radius_and_avg_height(l, h_avg);
    // If l is zero, then average height is simply min height
    if (std::fabs(l) < 1.0e-8)
    {
      h_avg = h_min;
    }

    // Loop over elements, disable refinement of elements above critical line
    unsigned n_el = Bulk_mesh_pt->nelement();
    for (unsigned e = 0; e < n_el; e++)
    {
      // Cast to element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

      // Find position of centre
      Vector<double> s(2, 0.0);
      Vector<double> x(2);
      el_pt->interpolated_x(s, x);

      // If outside a certain shell, disable refinement
      double r =
        pow((x[0] - l) * (x[0] - l) + (x[1] - h_avg) * (x[1] - h_avg), 0.5);
      double r1 = 0.3; // 0.4;
      double r2 = 0.2; // 0.3;
      double r3 = 0.1; // 0.2;
      double r4 = 0.05; // 0.1;
      double r5 = 0.025; // 0.05;
      double r6 = 0.0125;

      if (r > r1)
      {
        el_pt->disable_refinement();
      }
      else if (r > r2 && el_pt->refinement_level() >=
                           GlobalSimSettings::min_refinement_level + 1)
      {
        el_pt->disable_refinement();
      }
      else if (r > r3 && el_pt->refinement_level() >=
                           GlobalSimSettings::min_refinement_level + 2)
      {
        el_pt->disable_refinement();
      }
      else if (r > r4 && el_pt->refinement_level() >=
                           GlobalSimSettings::min_refinement_level + 3)
      {
        el_pt->disable_refinement();
      }
      else if (r > r5 && el_pt->refinement_level() >=
                           GlobalSimSettings::min_refinement_level + 4)
      {
        el_pt->disable_refinement();
      }
      else if (r > r6 && el_pt->refinement_level() >=
                           GlobalSimSettings::min_refinement_level + 5)
      {
        el_pt->disable_refinement();
      }
      else if (el_pt->refinement_level() >=
               GlobalSimSettings::min_refinement_level + 6)
      {
        el_pt->disable_refinement();
      }
      else
      { // Otherwise make sure its enabled
        el_pt->enable_refinement();
      }
    }
  }

  // Get the minimum height of the drop
  double get_min_height()
  {
    double h_min = 10000.0;
    unsigned n_nod = Bulk_mesh_pt->nboundary_node(1);
    for (unsigned n = 0; n < n_nod; n++)
    {
      double h_tmp = Bulk_mesh_pt->boundary_node_pt(1, n)->x(1);
      if (h_tmp < h_min)
      {
        h_min = h_tmp;
      }
    }
    std::cout << "h_min is: " << h_min << std::endl;
    return h_min;
  }

  /// Remap coordinates after adapting
  void remap_coordinates(Mesh* tmp_mesh_pt)
  {
    unsigned n_node = tmp_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      Node* nod_pt = tmp_mesh_pt->node_pt(n);

      // First move to origin
      nod_pt->x(0) -= GlobalSimSettings::h_initial_sim;

      // Rotate and move back
      double r_tmp = nod_pt->x(0);
      nod_pt->x(0) = nod_pt->x(1);
      nod_pt->x(1) = -r_tmp + GlobalSimSettings::h_initial_sim;
    }
  }

  /// Compute total energy of the elastic solid (gravitational, kinetic and
  /// elastic components)
  void compute_total_energy(double& E_g,
                            double& E_k,
                            double& E_e,
                            double& p_surf);

  /// Compute contact radius and average height in the film
  void compute_contact_radius_and_avg_height(double& l, double& h_avg);

  /// Timestepping routine
  void unsteady_run();

private:
  /// Create the boundary elements
  void create_boundary_elements();

  /// Delete boundary elements and wipe the associated mesh
  void delete_boundary_elements()
  {
    // How many surface elements are in the surface mesh
    unsigned n_element = Surface_mesh_pt->nelement();

    // Loop over the surface elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Kill surface element
      delete Surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    Surface_mesh_pt->flush_element_and_node_storage();
  } // end of delete_boundary_elements

  /// Complete the problem setup (attach function/physical parameter pointers)
  void complete_problem_setup();

  /// Bulk mesh
  ElasticRefineableHalfCircleSectorMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;

  /// Pointer to constitutive law
  ConstitutiveLaw* Constitutive_law_pt;

  // Pointers for mesh construction
  Circle* outer_boundary_circle_pt;

  // Pointer to error estimator
  Z2ErrorEstimator* error_estimator_pt;

  // DocInfo object
  DocInfo doc_info;
};


//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  ElasticLeidenfrostProblem()
{
#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  linear_solver_pt() = new MumpsSolver;
#endif

  // Set output directory
  doc_info.set_directory(GlobalSimSettings::result_folder);

  // Initialise counter for solutions
  doc_info.number() = 0;

  // Allocate the timestepper (this constructs the time object as well)
  // Boolean flag allows for temporal adaptivity
  add_time_stepper_pt(new TIMESTEPPER());

  // Don't reject timesteps above tolerance
  Keep_temporal_error_below_tolerance = false;

  // Circle object used to define the outer boundary
  outer_boundary_circle_pt = new Circle(0.0, 0.0, 1.0);

  // Create the mesh
  double x_centr = GlobalSimSettings::h_initial_sim;
  Bulk_mesh_pt = new ElasticRefineableHalfCircleSectorMesh<ELEMENT>(
    outer_boundary_circle_pt, x_centr, this->time_stepper_pt());

  // Set error estimator for bulk mesh
  error_estimator_pt = new Z2ErrorEstimator;

  // Set targets
  Bulk_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;
  Bulk_mesh_pt->max_permitted_error() = GlobalSimSettings::max_permitted_error;
  Bulk_mesh_pt->min_permitted_error() = GlobalSimSettings::min_permitted_error;
  Bulk_mesh_pt->max_refinement_level() =
    GlobalSimSettings::max_refinement_level;
  Bulk_mesh_pt->min_refinement_level() =
    GlobalSimSettings::min_refinement_level;

  // Create the "surface mesh" that will contain only the interface
  // elements. The constructor just creates the mesh without giving
  // it any elements, nodes, etc.
  Surface_mesh_pt = new Mesh;

  // Create the boundary elements
  create_boundary_elements();

  // Add the two submeshes
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all sub-meshes into a single mesh
  build_global_mesh();

  // Define a constitutive law: generalised Hookean
  // Constitutive_law_pt = new GeneralisedHookean(&GlobalPhysicalVariables::nu);

  StrainEnergyFunction* Strain_energy_function_pt;
  if (GlobalSimSettings::incompressible)
  {
    Strain_energy_function_pt =
      new IncompressibleNeoHookean(&GlobalPhysicalVariables::nu);
  }
  else
  {
    Strain_energy_function_pt = new NeoHookean(&GlobalPhysicalVariables::nu);
  }
  Constitutive_law_pt =
    new IsotropicStrainEnergyFunctionConstitutiveLaw(Strain_energy_function_pt);

  // Complete the problem set up
  complete_problem_setup();

  // Refine to get good resolution on the lower boundary
  for (unsigned n = 0; n < 5; n++)
  {
    adapt();
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      SolidNode* nod_pt = Bulk_mesh_pt->node_pt(n);
      nod_pt->x(0) = nod_pt->lagrangian_position(0);
      nod_pt->x(1) = nod_pt->lagrangian_position(1);
    }
  }

  // Now remap the coords, set the flag to true and do one final adapt to set up
  // vapour pressure elements correctly
  remap_coordinates(Bulk_mesh_pt);
  Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  GlobalSimSettings::base_mesh_created = true;
  adapt();

} // end constructor

//========start_of_complete_problem_setup==================================
/// Assign function/parameter pointers, set boundary conditions etc.
//=========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  complete_problem_setup()
{
  //----------------------Set up pointers in bulk mesh-----------------------//
  // Determine number of bulk elements in mesh
  const unsigned n_element_bulk = Bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e = 0; e < n_element_bulk; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the lambda squared parameter
    el_pt->lambda_sq_pt() = &GlobalPhysicalVariables::lambda_sq;

    // Set the forcing function
    el_pt->body_force_fct_pt() = &GlobalPhysicalVariables::get_body_force;

    // Set the constitutive law
    el_pt->constitutive_law_pt() = Constitutive_law_pt;

    // Set the element to be incompressible/compressible
    // Does nothing if not using a pressure formulation
    if (GlobalSimSettings::incompressible)
    {
      el_pt->set_incompressible();
    }
    else
    {
      el_pt->set_compressible();
    }

  } // End of loop over bulk elements
  //-------------------------------------------------------------------------//

  //----------------------Set up pointers in surface mesh--------------------//

  // Determine number of 1D interface elements in mesh
  const unsigned n_interface_element = Surface_mesh_pt->nelement();

  // Loop over the interface elements
  for (unsigned e = 0; e < n_interface_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    INTERFACE_ELEMENT* el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(e));

    // Set the pressure scaling parameter
    el_pt->pressure_scaling_param_pt() =
      &GlobalPhysicalVariables::pressure_scaling_param;

    // Set the conduction source term
    el_pt->conduction_source_pt() = &GlobalPhysicalVariables::conduction_source;

    // Set the film cutoff
    el_pt->film_cutoff_pt() = &GlobalSimSettings::film_cutoff;

    // Set the artificial repulsion
    el_pt->artificial_repulsion_fct_pt() =
      &GlobalPhysicalVariables::get_artificial_repulsion;

    // Get the outer unit normal - this is hacky, copying the output function
    Vector<double> normal(2);
    Vector<double> s(1);
    // Get local coordinates of plot point
    el_pt->get_s_plot(1, 3, s);
    // Get the outer unit normal
    el_pt->outer_unit_normal(s, normal);

    // Is this element part of the film?
    bool film_flag = false;
    if (normal[1] < -GlobalSimSettings::film_cutoff)
    {
      film_flag = true;
    }
    // Pin the gas pressure if not in the film
    const unsigned n_node = el_pt->nnode();
    for (unsigned n = 0; n < n_node; n++)
    {
      el_pt->unpin_vapour_pressure(n);
      if (film_flag == false)
      {
        el_pt->fix_vapour_pressure(n, 0.0);
      }
    }
  } // End of loop over interface elements

  //-------------------------------------------------------------------------//

  // Set the boundary conditions for this problem

  // Find number of nodes on boundary
  const unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(0);
  // Loop over nodes on fixed boundary
  for (unsigned n = 0; n < n_boundary_node; n++)
  {
    // Pin radial position
    Bulk_mesh_pt->boundary_node_pt(0, n)->pin_position(0);
  }
}

//============start_of_create_boundary_elements===============
/// Create boundary elements
//=======================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  create_boundary_elements()
{
  unsigned n_element = Bulk_mesh_pt->nboundary_element(1);
  // Loop over those elements adjacent to the free surface
  for (unsigned e = 0; e < n_element; e++)
  {
    // Set a pointer to the bulk element we wish to our interface
    // element to
    ELEMENT* bulk_element_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(1, e));

    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(1, e);

    // Create the interface element
    INTERFACE_ELEMENT* interface_element_pt =
      new INTERFACE_ELEMENT(bulk_element_pt, face_index);

    // Add the interface element to the surface mesh
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    interface_element_pt->set_boundary_number_in_bulk_mesh(1);
  }
} // end of create_boundary_elements


//==start_of_set_initial_condition========================================
/// Set initial conditions: Set all nodal velocities to zero and assign
/// acceleration
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  set_initial_condition()
{
  // First initialise the dts and set the weights
  initialise_dt(GlobalSimSettings::max_dt);

  double t_0;
  // Has there been any gravity set?
  if (GlobalPhysicalVariables::g_non_dim == 0.0)
  {
    t_0 =
      (GlobalPhysicalVariables::h_initial - GlobalSimSettings::h_initial_sim) /
      GlobalPhysicalVariables::W_non_dim;
  }
  else
  {
    t_0 = (-1.0 * GlobalPhysicalVariables::W_non_dim +
           pow(GlobalPhysicalVariables::W_non_dim *
                   GlobalPhysicalVariables::W_non_dim +
                 2.0 * GlobalPhysicalVariables::g_non_dim *
                   (GlobalPhysicalVariables::h_initial -
                    GlobalSimSettings::h_initial_sim),
               0.5)) /
          GlobalPhysicalVariables::g_non_dim;
  }

  // Set the initial condition manually
  // Number of steps in Newmark scheme
  const unsigned nsteps = 1;
  unsigned n_nod = Bulk_mesh_pt->nnode();
  for (unsigned n = 0; n < n_nod; n++)
  {
    // Cast to node
    Node* nod_pt = dynamic_cast<Node*>(Bulk_mesh_pt->node_pt(n));
    // Set previous nodal positions
    for (unsigned t = 1; t <= nsteps; t++)
    {
      double dt = t * GlobalSimSettings::max_dt;
      nod_pt->x(t, 0) = nod_pt->x(0, 0);
      nod_pt->x(t, 1) =
        nod_pt->x(0, 1) + GlobalPhysicalVariables::W_non_dim * dt +
        GlobalPhysicalVariables::g_non_dim * (t_0 * dt - 0.5 * dt * dt);
    }
    // Set the velocity and acceleration at the previous timestep
    nod_pt->x(nsteps + 1, 1) = -1.0 * (GlobalPhysicalVariables::W_non_dim +
                                       GlobalPhysicalVariables::g_non_dim *
                                         (t_0 - GlobalSimSettings::max_dt));
    nod_pt->x(nsteps + 2, 1) = -GlobalPhysicalVariables::g_non_dim;
  }
} // End of set_initial_condition


//===============start_of_custom_adaptive_timestepper=====================
/// A slight variation of adaptive_unsteady_newton_solve. Does a doc when
/// the solver fails
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  custom_adaptive_unsteady_newton_solve(const double& dt_desired,
                                        const double& epsilon,
                                        double& dt_new,
                                        bool& success)
{
  // Assume we will succeed
  success = true;

  // First, we need to backup the existing dofs, in case the timestep is
  // rejected

  // Find total number of dofs on current processor
  unsigned n_dof_local = dof_distribution_pt()->nrow_local();

  // Now set up a Vector to hold current values
  Vector<double> dofs_current(n_dof_local);

  // Load values into dofs_current
  for (unsigned i = 0; i < n_dof_local; i++) dofs_current[i] = dof(i);

  // Store the time
  double time_current = time_pt()->time();

  // Flag to detect whether the timestep has been rejected or not
  bool reject_timestep = 0;

  // The value of the actual timestep, by default the same as desired timestep
  double dt_actual = dt_desired;

  // Shift the time_values
  shift_time_values();

  // This loop surrounds the adaptive time-stepping and will not be broken
  // until a timestep is accepted
  do
  {
    // Initially we assume that this step will succeed and that this dt
    // value is ok.
    reject_timestep = 0;
    double dt_rescaling_factor = 1.0;

    // Set the new time and value of dt
    time_pt()->time() += dt_actual;
    time_pt()->dt() = dt_actual;

    // Set the predictor weights for the time stepper
    time_stepper_pt()->set_weights();

    // Run the timesteppers actions before timestep. These need to
    // be before the problem's actions_before_implicit_timestep so that the
    // boundary conditions are set consistently.
    time_stepper_pt()->actions_before_timestep(this);

    // Do any updates/boundary conditions changes here
    actions_before_implicit_timestep();

    // Attempt to solve the non-linear system
    try
    {
      // Solve the non-linear problem at this timestep
      newton_solve();
    }
    // Catch any exceptions thrown
    catch (NewtonSolverError& error)
    {
      // If it's a solver error then die
      if (error.linear_solver_error || Time_adaptive_newton_crash_on_solve_fail)
      {
        std::string error_message = "USER-DEFINED ERROR IN NEWTON SOLVER\n";
        error_message += "ERROR IN THE LINEAR SOLVER\n";

        // Die
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Reject the timestep, if we have an exception
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;

        // Half the time step
        dt_rescaling_factor = Timestep_reduction_factor_after_nonconvergence;
      }
    }

    // Run the timesteppers actions after timestep
    time_stepper_pt()->actions_after_timestep(this);

    // Update anything that needs updating after the timestep
    actions_after_implicit_timestep();

    // If timestep not already rejected, get an estimate for the overall vapour
    // pressure change. If this change is too large, reject the timestep
    if (!reject_timestep)
    {
      double p_avg = 0.0;
      double p_avg_m = 0.0;
      unsigned n_nod = Bulk_mesh_pt->nboundary_node(1);
      for (unsigned n = 0; n < n_nod; n++)
      {
        p_avg += Bulk_mesh_pt->boundary_node_pt(1, n)->value(0) / n_nod;
        p_avg_m += Bulk_mesh_pt->boundary_node_pt(1, n)->value(1, 0) / n_nod;
      }
      double p_var = std::fabs((p_avg - p_avg_m) / p_avg);
      std::cout << "Presssure variation: " << p_var << std::endl;
      if (p_var > GlobalSimSettings::epsilon_p)
      {
        oomph_info << "Vapour pressure variation: " << p_var
                   << " exceeds tolerance: " << GlobalSimSettings::epsilon_p
                   << std::endl;
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;
        dt_rescaling_factor /= 2.0;
      }
    }

    // If we haven't already rejected the timestep, do a check for inverted
    // elements
    if (!reject_timestep)
    {
      bool has_inverted_elements = false;
      Bulk_mesh_pt->check_inverted_elements(has_inverted_elements);
      if (has_inverted_elements)
      {
        // Set the flag and half the new timestep candidate
        oomph_info << "Inverted element(s) detected: " << std::endl;
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;
        dt_rescaling_factor /= 2.0;
      }
    }

    // If we're not rejecting the timestep
    // then calculate the error estimate and rescaling factor.
    if (!reject_timestep)
    {
      // Get a global error norm to use in adaptivity (as specified by the
      // problem sub-class writer). Prevent a divide by zero if the solution
      // gives very close to zero error. Error norm should never be negative
      // but use absolute value just in case.
      double error =
        std::max(std::abs(custom_global_temporal_error_norm()), 1e-12);

      // Calculate the scaling  factor
      dt_rescaling_factor = pow(epsilon / error, 0.5);

      oomph_info << "Timestep scaling factor is  " << dt_rescaling_factor
                 << std::endl;
      oomph_info << "Estimated timestepping error is " << error << std::endl;

      // Do we have to do it again?
      if (error > epsilon)
      {
        oomph_info << "Estimated timestepping error " << error
                   << " exceeds tolerance " << epsilon << " ";
        if (Keep_temporal_error_below_tolerance)
        {
          oomph_info << " --> rejecting timestep." << std::endl;
          reject_timestep = 1;
        }
        else
        {
          oomph_info << " ...but we're not rejecting the timestep" << std::endl;
        }
        oomph_info
          << "Note: This behaviour can be adjusted by changing the protected "
          << "boolean" << std::endl
          << std::endl
          << "    Problem::Keep_temporal_error_below_tolerance" << std::endl;
      }
    } // End of !reject_timestep flag

    // Calculate the next time step size and check it's ok
    // ============================================================

    // Calculate the possible next time step, if no error conditions
    // trigger.
    double new_dt_candidate = dt_rescaling_factor * dt_actual;
    double custom_time_scaling_factor = 2.0;

    // Check that the scaling factor is within the allowed range
    if (dt_rescaling_factor > custom_time_scaling_factor)
    {
      oomph_info << "Tried to increase dt by the ratio " << dt_rescaling_factor
                 << " which is above the maximum ("
                 << custom_time_scaling_factor
                 << "). Attempting to increase by the maximum ratio instead."
                 << std::endl;
      new_dt_candidate = custom_time_scaling_factor * dt_actual;
    }
    // If we have already rejected the timestep then don't do this check
    // because DTSF will definitely be too small.
    else if ((!reject_timestep) && (dt_rescaling_factor <= DTSF_min_decrease))
    {
      // Handle this special case where we want to continue anyway (usually
      // Minimum_dt_but_still_proceed = -1 so this has no effect).
      if (new_dt_candidate < Minimum_dt_but_still_proceed)
      {
        oomph_info
          << "Warning: Adaptation of timestep to ensure satisfaction \n "
          << "of error bounds during adaptive timestepping \n "
          << "would lower dt below \n"
          << "Problem::Minimum_dt_but_still_proceed="
          << Minimum_dt_but_still_proceed << "\n"
          << "---> We're continuing with present timestep.\n " << std::endl;
        dt_rescaling_factor = 1.0;
      }
      else
      {
        // Otherwise reject
        oomph_info << "Timestep would decrease by " << dt_rescaling_factor
                   << " which is less than the minimum scaling factor "
                   << DTSF_min_decrease << std::endl;
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;
      }
    }

    // Now check that the new dt is within the allowed range
    if (new_dt_candidate > Maximum_dt)
    {
      oomph_info << "Tried to increase dt to " << new_dt_candidate
                 << " which is above the maximum (" << Maximum_dt
                 << "). I increased it to the maximum value instead.";
      dt_actual = Maximum_dt;
    }
    else if (new_dt_candidate < Minimum_dt)
    {
      oomph_info << "MINIMUM TIMESTEP HAS BEEN REACHED: BREAKING OUT AND DOING "
                    "A FINAL FULL DOC."
                 << std::endl;
      doc_solution();
      success = false;
      break;
    }
    else
    {
      dt_actual = new_dt_candidate;
    }

    actions_after_implicit_timestep_and_error_estimation();

    // If we are rejecting this attempt then revert the dofs etc.
    if (reject_timestep)
    {
      // Reset the time
      time_pt()->time() = time_current;

      // Reload the dofs
      unsigned ni = dofs_current.size();
      for (unsigned i = 0; i < ni; i++)
      {
        dof(i) = dofs_current[i];
      }

#ifdef OOMPH_HAS_MPI
      // Synchronise the solution on different processors (on each submesh)
      this->synchronise_all_dofs();
#endif

      // Call all "after" actions, e.g. to handle mesh updates
      actions_after_newton_step();
      actions_before_newton_convergence_check();
      actions_after_newton_solve();
      actions_after_implicit_timestep();
      actions_after_implicit_timestep_and_error_estimation();
    }

  }
  // Keep this loop going until we accept the timestep
  while (reject_timestep);

  // Once the timestep has been accepted, return the time step that should be
  // used next time.
  dt_new = dt_actual;
}

//===============start_of_custom_global_temporal_error_norm===============
/// Temporal error estimation for the Newmark scheme. Uses a more accurate
/// interpolation of the acceleration between timesteps
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
double ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  custom_global_temporal_error_norm()
{
  // First calculate the centre of mass
  unsigned n_el = Bulk_mesh_pt->nelement();
  double z_com_current = 0.0;
  double vol = 0.0;
  for (unsigned e = 0; e < n_el; e++)
  {
    // Cast to element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Loop over integration points
    unsigned n_intpt = el_pt->integral_pt()->nweight();
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Jacobian of mapping
      double J = el_pt->J_eulerian_at_knot(ipt);

      // Integral weight
      double w = el_pt->integral_pt()->weight(ipt);

      // Geometric Jacobian
      Vector<double> s(2);
      for (unsigned i = 0; i < 2; ++i)
      {
        s[i] = el_pt->integral_pt()->knot(ipt, i);
      }
      Vector<double> x(2, 0.0);
      el_pt->interpolated_x(s, x);

      // Add contribution
      z_com_current += x[1] * x[0] * J * w;
      vol += x[0] * J * w;
    }
  }
  z_com_current /= vol;

  // Get the timestep from the last solve
  double dt_local = time_pt()->dt();

  // The beta2 parameter from the Newmark timestepper (probably 0.5)
  double beta = GlobalSimSettings::Newmark_beta2;

  // Number of nodes in the mesh
  unsigned nnod = Bulk_mesh_pt->nnode();

  // Allocate average error
  double relative_error_avg = 0.0;

  // Number of steps used for the Newmark timestepper
  const unsigned nsteps = 1;

  // Do a first loop over the nodes to find the largest radial and vertical
  // displacements
  double disp_r_max = 0.0;
  double disp_z_max = 0.0;
  for (unsigned n = 0; n < nnod; n++)
  {
    // Current displacements
    double disp_r = std::fabs(Bulk_mesh_pt->node_pt(n)->position(0) -
                              Bulk_mesh_pt->node_pt(n)->lagrangian_position(0));
    double disp_z =
      std::fabs(Bulk_mesh_pt->node_pt(n)->position(1) -
                Bulk_mesh_pt->node_pt(n)->lagrangian_position(1) -
                (z_com_current - GlobalSimSettings::h_initial_sim));

    if (disp_r > disp_r_max)
    {
      disp_r_max = disp_r;
    }
    if (disp_z > disp_z_max)
    {
      disp_z_max = disp_z;
    }
  }

  // Loop over nodes and estimate error
  for (unsigned n = 0; n < nnod; n++)
  {
    // Cast to the node
    Node* nod_pt = dynamic_cast<Node*>(Bulk_mesh_pt->node_pt(n));

    // The accelerations at previous time
    double r_ddot_prev = nod_pt->position(nsteps + 2, 0);
    double z_ddot_prev = nod_pt->position(nsteps + 2, 1);

    // The current accelerations
    double r_ddot_current = nod_pt->dposition_dt(2, 0);
    double z_ddot_current = nod_pt->dposition_dt(2, 1);

    // The error components
    double abs_error_r = std::fabs((-1.0 / 6.0 + beta) * dt_local * dt_local *
                                   (r_ddot_current - r_ddot_prev));
    double abs_error_z = std::fabs((-1.0 / 6.0 + beta) * dt_local * dt_local *
                                   (z_ddot_current - z_ddot_prev));

    // Add to the relative error average
    relative_error_avg +=
      0.5 * (abs_error_r / disp_r_max + abs_error_z / disp_z_max);
  }

  // Divide by the number of nodes
  relative_error_avg /= nnod;

  // Now look at interface contribution, get average height and contact radius
  double l, h_avg;
  compute_contact_radius_and_avg_height(l, h_avg);

  // Return the error
  return relative_error_avg;
}

//===============start_of_compute_total_energy============================
/// Compute the total energy for the elastic solid (gravitational +
/// kinetic + elastic)
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  compute_total_energy(double& E_g, double& E_k, double& E_e, double& p_surf)
{
  // Number of elements in the bulk mesh;
  unsigned n_el = Bulk_mesh_pt->nelement();

  // Initialize values
  E_g = 0.0;
  E_k = 0.0;
  E_e = 0.0;

  // Container for local coordinates
  Vector<double> s(2);

  // Loop over elements
  for (unsigned e = 0; e < n_el; e++)
  {
    // Cast to element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Number of nodes in the element
    unsigned n_node = el_pt->nnode();

    // Container for shape function
    Shape psi(n_node);
    DShape dpsidxi(n_node, 2);

    // Find the number integration points
    unsigned n_intpt = el_pt->integral_pt()->nweight();

    // Loop over intregration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign the values of s
      for (unsigned i = 0; i < 2; ++i)
      {
        s[i] = el_pt->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = el_pt->integral_pt()->weight(ipt);

      // Call the shape funcitons and get the Jacobian
      double J = el_pt->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

      // Containers for Lagrangian and Eulerian position
      Vector<double> interpolated_xi(2, 0.0);
      Vector<double> interpolated_X(2, 0.0);
      DenseMatrix<double> interpolated_dXdxi(2, 2, 0.0);
      Vector<double> veloc(2, 0.0);

      // Loop over nodes to calculate displacement and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the nodal displacements
        for (unsigned i = 0; i < 2; i++)
        {
          interpolated_xi[i] += el_pt->lagrangian_position(l, i) * psi(l);
          interpolated_X[i] += el_pt->nodal_position(l, i) * psi(l);

          // Loop over Lagrangian derivative directions
          for (unsigned j = 0; j < 2; j++)
          {
            // Calculate dX[i]/dxi_{j}
            interpolated_dXdxi(i, j) +=
              el_pt->nodal_position(l, i) * dpsidxi(l, j);
          }

          // Calculate the velocity components
          veloc[i] += el_pt->dnodal_position_dt(l, i) * psi(l);
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

      // Initialise
      DenseMatrix<double> gup(3, 3, 0.0);
      DenseMatrix<double> Gup(3, 3, 0.0);
      gup(0, 0) = 1.0;
      gup(1, 1) = 1.0;
      gup(2, 2) = 1.0 / (interpolated_xi[0] * interpolated_xi[0]);
      double detg = interpolated_xi[0] * interpolated_xi[0];
      double detG = interpolated_X[0] * interpolated_X[0] *
                    ((interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 0) +
                      interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 0)) *
                       (interpolated_dXdxi(0, 1) * interpolated_dXdxi(0, 1) +
                        interpolated_dXdxi(1, 1) * interpolated_dXdxi(1, 1)) -
                     pow(interpolated_dXdxi(0, 0) * interpolated_dXdxi(0, 1) +
                           interpolated_dXdxi(1, 0) * interpolated_dXdxi(1, 1),
                         2.0));
      Vector<double> I(3, 0.0);
      I[2] = detG / detg;
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          I[0] += gup(i, j) * G(i, j);
          I[1] += Gup(i, j) * g(i, j) * I[2];
        }
      }

      // Pre-multiply weights, jacobian and the square root of the determinant
      // of the undeformed metric tensor (xi[0])
      double W = w * J * interpolated_xi[0];

      // Add integral contributions
      // Gravity
      E_g += GlobalPhysicalVariables::g_non_dim * interpolated_X[1] * W;
      // Elastic
      E_e += 0.5 *
             (0.5 * (I[0] - 3.0) - 0.5 * log(I[2]) +
              (GlobalPhysicalVariables::nu /
               (1.0 - 2.0 * GlobalPhysicalVariables::nu)) *
                0.25 * log(I[2]) * log(I[2])) /
             (1.0 + GlobalPhysicalVariables::nu) * W;
      // Kinetic
      E_k += GlobalPhysicalVariables::lambda_sq * 0.5 *
             (veloc[0] * veloc[0] + veloc[1] * veloc[1]) * W;
    }
  }

  // Now loop over surface elements, to calculate the power from external
  // forcing
  p_surf = 0.0;
  unsigned n_surface_el = Surface_mesh_pt->nelement();
  for (unsigned e = 0; e < n_surface_el; e++)
  {
    // Cast to element
    INTERFACE_ELEMENT* el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(e));

    double p_el;
    el_pt->compute_surface_power_contribution(p_el);
    p_surf += p_el;
  }
}

//===============start_of_compute_contact_radius_and_avg_height===========
/// Compute the contact radius (during impact) and average height by
/// integrating up to the contact radius
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  compute_contact_radius_and_avg_height(double& l, double& h_avg)
{
  // // Initialize integrals
  // double h_int = 0.0;
  // double area_int = 0.0;

  // // Number of elements in the surface mesh;
  // unsigned n_el = Surface_mesh_pt->nelement();

  // // Boolean for inclusion of integral contributions
  // bool first = true;

  // // Loop over surface elements
  // for (int e = n_el - 1; e >= 0; e--)
  // {
  //   // Cast to element
  //   INTERFACE_ELEMENT* el_pt =
  //     dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(e));

  //   // Get the outer unit normal - this is hacky, copying the output function
  //   Vector<double> normal(2);
  //   Vector<double> x(2, 0.0);
  //   Vector<double> s(1);
  //   // Get local coordinates of plot point
  //   el_pt->get_s_plot(1, 3, s);
  //   // Get the position vector
  //   el_pt->interpolated_x(s, x);
  //   // Get the outer unit normal
  //   el_pt->outer_unit_normal(s, normal);

  //   // First time the x component of the normal changes sign: change the
  //   // Boolean
  //   if (normal[0] < 0.0)
  //   {
  //     first = false;
  //   }

  //   // If the element is in the film add the integral contribution
  //   if (!first)
  //   {
  //     // Create a new QElement
  //     QElement<1, 3>* Qel_pt = new QElement<1, 3>();

  //     // Set the nodal positions of the QElement (radial nodal position of
  //     the
  //     // surface element)
  //     for (unsigned l = 0; l < 3; l++)
  //     {
  //       Qel_pt->node_pt(l) = new Node(1, 1, 1, true);
  //       Qel_pt->node_pt(l)->x(0) = el_pt->nodal_position(l, 0);
  //     }

  //     // Number of integration points
  //     unsigned n_intpt = Qel_pt->integral_pt()->nweight();

  //     // Set up memory for shape functions
  //     Shape psi_v(3);
  //     DShape dpsi_vdx(3, 1);

  //     // Loop over integration points
  //     for (unsigned ipt = 0; ipt < n_intpt; ipt++)
  //     {
  //       // Get the integral weight
  //       double w = Qel_pt->integral_pt()->weight(ipt);

  //       // Call the shape functions and derivatives
  //       double J =
  //         std::fabs(Qel_pt->dshape_eulerian_at_knot(ipt, psi_v, dpsi_vdx));

  //       // Positions
  //       Vector<double> x2(2, 0.0);
  //       for (unsigned l = 0; l < 3; l++)
  //       {
  //         for (unsigned i = 0; i < 2; i++)
  //         {
  //           // Positions
  //           x2[i] += el_pt->nodal_position(l, i) * psi_v(l);
  //         }
  //       }

  //       // Add integral contributions
  //       h_int += x2[0] * x2[1] * w * J;
  //       area_int += x2[0] * w * J;
  //     }

  //     // Delete allocated QElement and its nodes
  //     for (unsigned l = 0; l < 3; l++)
  //     {
  //       delete Qel_pt->node_pt(l);
  //     }
  //     delete Qel_pt;
  //   }
  // }

  // // Calculate values
  // l = pow(2.0 * area_int, 0.5);
  // h_avg = h_int / area_int;

  // Number of interior points per surface element (at which the normal is
  // sampled)
  unsigned n_interior = 4;

  // Tolerance for variation in the normal
  double eps_normal = 0.025;

  // Number of surface elements
  unsigned n_el = Surface_mesh_pt->nelement();

  // Initialize height integral
  double h_int = 0.0;

  // Previous position
  Vector<double> x_prev(2, 0.0);

  DenseMatrix<double> r_vals(n_el, n_interior, 0.0);
  DenseMatrix<double> r_vals_prev(n_el, n_interior, 0.0);
  DenseMatrix<double> z_vals(n_el, n_interior, 0.0);
  DenseMatrix<double> z_vals_prev(n_el, n_interior, 0.0);

  // Loop over surface elements
  for (unsigned e = 0; e < n_el; e++)
  {
    // Cast to element
    INTERFACE_ELEMENT* el_pt =
      dynamic_cast<INTERFACE_ELEMENT*>(Surface_mesh_pt->element_pt(e));

    // Move the query points off the element boundaries
    bool use_equally_spaced_interior_sample_points = true;

    // Containers for normal etc.
    Vector<double> normal(2);
    Vector<double> x(2, 0.0);
    Vector<double> s(1);

    // Loop over the interior points
    for (unsigned n = 0; n < n_interior; n++)
    {
      // Get the outer unit normal - this is hacky, copying the output function
      // Get local coordinates of plot point
      el_pt->get_s_plot(
        n, n_interior, s, use_equally_spaced_interior_sample_points);
      // Get the position vector
      el_pt->interpolated_x(s, x);
      // Get the outer unit normal
      el_pt->outer_unit_normal(s, normal);

      // Assign coords
      r_vals(e, n) = x[0];
      r_vals_prev(e, n) = x_prev[0];
      z_vals(e, n) = x[1];
      z_vals_prev(e, n) = x_prev[1];

      // Assign the neck value if within tolerance
      double diff =
        pow(normal[0] * normal[0] + (normal[1] + 1.0) * (normal[1] + 1.0), 0.5);
      if (diff < eps_normal && x[0] < 0.6)
      {
        l = x[0];
      }

      // Re-assign previous position
      x_prev = x;
    }
  }

  // Loop again to calculate height integral
  for (unsigned e = 0; e < n_el; e++)
  {
    for (unsigned n = 0; n < n_interior; n++)
    {
      if (r_vals(e, n) <= l && z_vals(e, n) < 0.5)
      {
        h_int +=
          (r_vals(e, n) - r_vals_prev(e, n)) * 0.5 *
          (r_vals(e, n) * z_vals(e, n) + r_vals_prev(e, n) * z_vals_prev(e, n));
      }
    }
  }

  // Calculate average height
  h_avg = (2.0 / (l * l)) * h_int;
}

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 3;

  if (doc_info.number() % GlobalSimSettings::steps_per_interface_doc == 0)
  {
    sprintf(filename,
            "%s/interface_soln%i.dat",
            doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    Surface_mesh_pt->output(some_file, npts);
    some_file.close();
  }

  // npts = 2;
  // if (doc_info.number() % GlobalSimSettings::steps_per_bulk_doc == 0)
  // {
  //   sprintf(filename,
  //           "%s/bulk_soln%i.dat",
  //           doc_info.directory().c_str(),
  //           doc_info.number());
  //   some_file.open(filename);
  //   Bulk_mesh_pt->output(some_file, npts);
  //   some_file.close();

  //   // Also do a paraview output
  //   // sprintf(filename,
  //   //         "%s/bulk_soln%i.vtu",
  //   //         doc_info.directory().c_str(),
  //   //         doc_info.number());
  //   // some_file.open(filename);
  //   // Bulk_mesh_pt->output_paraview(some_file, npts);
  //   // some_file.close();
  // }

  // Get the total energy
  double E_g, E_k, E_e, p_surf;
  compute_total_energy(E_g, E_k, E_e, p_surf);

  // Get the contact radius and avg height
  double l, h_avg;
  compute_contact_radius_and_avg_height(l, h_avg);

  // Current time
  ofstream timefile;
  timefile.open(GlobalSimSettings::result_folder + "/time.dat",
                std::ofstream::out | std::ofstream::app);
  timefile << time_pt()->time() << " " << E_g << " " << E_k << " " << E_e << " "
           << p_surf << " " << l << " " << h_avg << "\n";
  timefile.close();

  // Increment the counter
  doc_info.number()++;
}


//========================================================================
/// Unsteady run
//========================================================================
template<class ELEMENT, class INTERFACE_ELEMENT, class TIMESTEPPER>
void ElasticLeidenfrostProblem<ELEMENT, INTERFACE_ELEMENT, TIMESTEPPER>::
  unsteady_run()
{
  // Start time of the run
  time_t start_time = std::time(NULL);

  set_initial_condition();
  doc_solution();

  // Initial timestep: Use the one used when setting up the initial
  // condition
  double dt = time_pt()->dt();

  double dt_new = dt;
  bool success = true;
  bool below_impact_height = false;
  unsigned n = 0;
  while (success)
  {
    // If we're running over time (47:00:00), do one final full doc and break
    // out of the loop
    time_t current_time = std::time(NULL);
    if (current_time - start_time > 21600.0)
    {
      oomph_info << "TIME LIMIT HAS BEEN REACHED. BREAKING OUT." << std::endl;
      break;
    }

    // Perform solve
    if (GlobalSimSettings::use_adaptive_timestepping)
    {
      custom_adaptive_unsteady_newton_solve(
        dt, GlobalSimSettings::epsilon_t, dt_new, success);
    }
    else
    {
      (void)adaptive_unsteady_newton_solve(dt, true);
    }

    // Set the timestep
    dt = dt_new;

    // Doc interface/bulk
    doc_solution();

    // Check whether temporal/spatial adaptivity should be used
    double h0 = get_min_height();
    toggle_bulk_mesh_refinement(h0);

    // Adapt if necessary
    if (n % GlobalSimSettings::steps_per_adapt == 0)
    {
      adapt();
      toggle_bulk_mesh_refinement(h0);
      adapt();
    }
    n++;

    if (h0 < 0.0025)
    {
      // Set the tolerance for pressure change
      GlobalSimSettings::epsilon_p = 3.0e-2;
    }

    // Check whether we need to stop the simulation (i.e. if rebound has been
    // detected)
    if (h0 < 0.01 && !below_impact_height)
    {
      // Switched to true when the drop first makes impact
      below_impact_height = true;
    }
    if (h0 > 0.01 && below_impact_height)
    {
      // The drop has just rebounded
      success = false;
      std::cout << "REBOUND DETECTED. STOPPING RUN." << std::endl;
    }
  }
}


//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc, argv);
  // oomph_mpi_output.restrict_output_to_single_processor();
#endif

  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--use_solid_pressure",
    &GlobalSimSettings::use_pressure_formulation,
    "Use solid pressure");

  CommandLineArgs::specify_command_line_flag(
    "--incompressible", &GlobalSimSettings::incompressible, "Incompressible");

  CommandLineArgs::specify_command_line_flag(
    "--folder", &GlobalSimSettings::result_folder, "Result Folder");

  CommandLineArgs::specify_command_line_flag(
    "--alpha",
    &GlobalPhysicalVariables::pressure_scaling_param,
    "Pressure scaling number");

  CommandLineArgs::specify_command_line_flag(
    "--poisson", &GlobalPhysicalVariables::nu, "Poisson ratio");

  CommandLineArgs::specify_command_line_flag(
    "--gravity", &GlobalPhysicalVariables::g_non_dim, "Gravity");

  CommandLineArgs::specify_command_line_flag(
    "--evap_number",
    &GlobalPhysicalVariables::conduction_source,
    "Evaporation number");

  CommandLineArgs::specify_command_line_flag(
    "--initial_height",
    &GlobalPhysicalVariables::h_initial,
    "Initial drop height");

  CommandLineArgs::specify_command_line_flag("--initial_height_sim",
                                             &GlobalSimSettings::h_initial_sim,
                                             "Simulation starting height");

  CommandLineArgs::specify_command_line_flag(
    "--initial_speed",
    &GlobalPhysicalVariables::W_non_dim,
    "Initial drop speed");

  CommandLineArgs::specify_command_line_flag(
    "--max_timestep", &GlobalSimSettings::max_dt, "Maximum timestep");

  CommandLineArgs::specify_command_line_flag("--newmark_beta_1",
                                             &GlobalSimSettings::Newmark_beta1,
                                             "First Newmark parameter");

  CommandLineArgs::specify_command_line_flag("--newmark_beta_2",
                                             &GlobalSimSettings::Newmark_beta2,
                                             "Second Newmark parameter");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // If the given start height is less than the sim start height, set them
  // equal
  if (GlobalPhysicalVariables::h_initial < GlobalSimSettings::h_initial_sim)
  {
    GlobalSimSettings::h_initial_sim = GlobalPhysicalVariables::h_initial;
  }

  // Clear timefile before start
  ofstream timefile(GlobalSimSettings::result_folder + "/time.dat");

  // Create the problem
  if (GlobalSimSettings::use_pressure_formulation)
  {
    // Build the problem with an additional degree of freedom for pressure
    ElasticLeidenfrostProblem<
      RefineableQAxisymmetricCylindricalPVDWithPressureElement,
      SolidLubricationElement<
        RefineableQAxisymmetricCylindricalPVDWithPressureElement>,
      VariableNewmark<1>>
      problem;

    problem.unsteady_run();
  }
  else
  {
    // Standard build
    ElasticLeidenfrostProblem<
      RefineableQAxisymmetricCylindricalPVDElement<3>,
      SolidLubricationElement<RefineableQAxisymmetricCylindricalPVDElement<3>>,
      VariableNewmark<1>>
      problem;

    problem.unsteady_run();
  }

  // Finalise MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end main
