// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#ifndef OOMPH_HALF_CIRCLE_SECTOR_MESH_TEMPLATE_CC
#define OOMPH_HALF_CIRCLE_SECTOR_MESH_TEMPLATE_CC


#include "half_circle_sector_mesh.template.h"

namespace oomph
{
  //====================================================================
  /// Constructor for deformable 2D Ring mesh class. Pass pointer to
  /// geometric object that specifies the wall, start and end coordinates on the
  /// geometric object, and the fraction along
  /// which the dividing line is to be placed, and the timestepper
  /// (defaults to (Steady) default timestepper defined in Mesh).
  /// Nodal positions are determined via macro-element-based representation
  /// of the Domain (as a QuarterCircleSectorDomain).
  //====================================================================
  template<class ELEMENT>
  HalfCircleSectorMesh<ELEMENT>::HalfCircleSectorMesh(
    GeomObject* wall_pt, double& x_centr, TimeStepper* time_stepper_pt)
    : Wall_pt(wall_pt), X_centr(x_centr)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Build macro element-based domain
    Domain_pt = new HalfCircleSectorDomain(wall_pt, x_centr);

    // Set the number of boundaries
    set_nboundary(2);

    // We have only bothered to parametrise boundary 1
    // Boundary_coordinate_exists[1] = true;

    // Allocate the store for the elements
    Element_pt.resize(4);

    // Create first element
    Element_pt[0] = new ELEMENT;

    // Read out the number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    // Can now allocate the store for the nodes
    Node_pt.resize(4 * n_p * n_p - 3 * n_p - 2 * (n_p - 1));


    Vector<double> s(2);
    Vector<double> r(2);

    // Storage for the intrinsic boundary coordinate
    Vector<double> zeta(1);


    // Set up geometrical data
    //------------------------

    // Initialise node counter
    unsigned long node_count = 0;


    // Now assign the topology
    // Boundaries are numbered 0 1 2 from the bottom proceeding anticlockwise


    // FIRST ELEMENT (centre)
    //

    // Set the corner node (on boundaries 0 and 2)
    //-------------------------------------------

    // Create the ll node
    Node_pt[node_count] =
      finite_element_pt(0)->construct_boundary_node(0, time_stepper_pt);

    // Set the pointer from the element to the node
    finite_element_pt(0)->node_pt(0) = Node_pt[node_count];

    // Set the position of the ll node
    s[0] = -1.0;
    s[1] = -1.0;
    Domain_pt->macro_element_pt(0)->macro_map(s, r);
    Node_pt[node_count]->x(0) = r[0];
    Node_pt[node_count]->x(1) = r[1];

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // First row is on boundary 0:
    //---------------------------
    for (unsigned l1 = 1; l1 < n_p; l1++)
    {
      // Local node number
      unsigned jnod_local = l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(0)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(0)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = -1.0;
      Domain_pt->macro_element_pt(0)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // Loop over the other rows of nodes
    //------------------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // The nodes are in the interior
      //------------------------------------
      // Loop over the other node columns
      for (unsigned l1 = 0; l1 < n_p; l1++)
      {
        // Local node number
        unsigned jnod_local = l1 + n_p * l2;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(0)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(0)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(0)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }

    // SECOND ELEMENT (lower right corner)
    // Create element
    Element_pt[1] = new ELEMENT;

    // Loop over the first column (already exists!)
    //---------------------------------------------
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      // Node number in existing element
      unsigned jnod_local_old = (n_p - 1) + l2 * n_p;

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(l2 * n_p) =
        finite_element_pt(0)->node_pt(jnod_local_old);
    }

    // Loop over the other node columns (apart from last one)
    //------------------------------------------------------
    for (unsigned l1 = 1; l1 < n_p - 1; l1++)
    {
      // First node is at the bottom (on boundary 0)
      //--------------------------------------------

      // Local node number
      unsigned jnod_local = l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = -1.0;
      Domain_pt->macro_element_pt(1)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now loop over the interior nodes in this column
      //-------------------------------------------------
      for (unsigned l2 = 1; l2 < n_p; l2++)
      {
        // Local node number
        unsigned jnod_local = l1 + l2 * n_p;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(1)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(1)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }

    // Last column (on boundary 1)
    //----------------------------

    // First node is at the bottom (and hence also on boundary 0)
    //-----------------------------------------------------------

    // Local node number
    unsigned jnod_local = n_p - 1;

    // Create the node
    Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
      jnod_local, time_stepper_pt);

    // Set the pointer from the element to the node
    finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

    // Set the position of the node
    s[0] = 1.0;
    s[1] = -1.0;
    Domain_pt->macro_element_pt(1)->macro_map(s, r);
    Node_pt[node_count]->x(0) = r[0];
    Node_pt[node_count]->x(1) = r[1];

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the remaining nodes in last column (only on boundary 1)
    //---------------------------------------------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Local node number
      unsigned jnod_local = (n_p - 1) + l2 * n_p;

      // Create the node
      Node_pt[node_count] = finite_element_pt(1)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(1)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = 1.0;
      s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
      Domain_pt->macro_element_pt(1)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }


    // THIRD ELEMENT (top)
    // Create element
    Element_pt[2] = new ELEMENT;

    // Loop over the first row (has already been created via element 0)
    //-----------------------------------------------------------------
    for (unsigned l1 = 0; l1 < n_p; l1++)
    {
      // Node number in existing element
      unsigned jnod_local_old = n_p * (n_p - 1) + l1;

      // Local node number here
      unsigned jnod_local = l1;

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) =
        finite_element_pt(0)->node_pt(jnod_local_old);
    }


    // Loop over the remaining nodes in the last column (has already
    //--------------------------------------------------------------
    // been created via element 1)
    //----------------------------
    for (unsigned l2 = 1; l2 < n_p; l2++)
    {
      // Node number in existing element
      unsigned jnod_local_old = n_p * (n_p - 1) + l2;

      // Local node number here
      unsigned jnod_local = (n_p - 1) + l2 * n_p;

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) =
        finite_element_pt(1)->node_pt(jnod_local_old);
    }


    // Loop over the nodes in rows (apart from last one which is on boundary 1)
    //-------------------------------------------------------------------------
    for (unsigned l2 = 1; l2 < n_p - 1; l2++)
    {
      // The nodes are in the interior
      //------------------------------------
      // Loop over the other node columns
      for (unsigned l1 = 0; l1 < n_p - 1; l1++)
      {
        // Local node number
        unsigned jnod_local = l1 + n_p * l2;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(2)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(2)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
    }

    // Top row is on boundary 1 only:
    //---------------------------------------
    for (unsigned l1 = 0; l1 < n_p - 1; l1++)
    {
      // Local node number
      unsigned jnod_local = n_p * (n_p - 1) + l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(2)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(2)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = 1.0;
      Domain_pt->macro_element_pt(2)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // FOURTH ELEMENT (bottom left)
    // Create element
    Element_pt[3] = new ELEMENT;

    // Loop over the last column (already exists!)
    //---------------------------------------------
    for (unsigned l2 = 0; l2 < n_p; l2++)
    {
      // Node number in existing element
      unsigned jnod_local_old = l2 * n_p;

      // Set the pointer from the element to the node
      finite_element_pt(3)->node_pt(l2 * n_p + (n_p - 1)) =
        finite_element_pt(0)->node_pt(jnod_local_old);
    }

    // Loop over the other node columns (apart from the first one)
    //------------------------------------------------------
    for (unsigned l1 = 1; l1 < n_p - 1; l1++)
    {
      // First node is at the bottom (on boundary 0)
      //--------------------------------------------

      // Local node number
      unsigned jnod_local = l1;

      // Create the node
      Node_pt[node_count] = finite_element_pt(3)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(3)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
      s[1] = -1.0;
      Domain_pt->macro_element_pt(3)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now loop over the interior nodes in this column
      //-------------------------------------------------
      for (unsigned l2 = 1; l2 < n_p - 1; l2++)
      {
        // Local node number
        unsigned jnod_local = l1 + l2 * n_p;

        // Create the node
        Node_pt[node_count] =
          finite_element_pt(3)->construct_node(jnod_local, time_stepper_pt);

        // Set the pointer from the element to the node
        finite_element_pt(3)->node_pt(jnod_local) = Node_pt[node_count];

        // Set the position of the node
        s[0] = -1.0 + 2.0 * double(l1) / double(n_p - 1);
        s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
        Domain_pt->macro_element_pt(3)->macro_map(s, r);
        Node_pt[node_count]->x(0) = r[0];
        Node_pt[node_count]->x(1) = r[1];

        // Increment the node number
        node_count++;
      }
      // Nodes on the top row have already been created and belong to the top
      // element

      // Local node number
      jnod_local = l1 + (n_p - 1) * n_p;

      // Node number in existing element
      unsigned jnod_local_old = n_p * l1;

      // Set the pointer from the element to the node
      finite_element_pt(3)->node_pt(jnod_local) =
        finite_element_pt(2)->node_pt(jnod_local_old);
    }
    for (unsigned l2 = 0; l2 < n_p - 1; l2++)
    {
      // Local node number
      unsigned jnod_local = n_p * l2;

      // Create the node
      Node_pt[node_count] = finite_element_pt(3)->construct_boundary_node(
        jnod_local, time_stepper_pt);

      // Set the pointer from the element to the node
      finite_element_pt(3)->node_pt(jnod_local) = Node_pt[node_count];

      // Set the position of the node
      s[0] = -1.0;
      s[1] = -1.0 + 2.0 * double(l2) / double(n_p - 1);
      Domain_pt->macro_element_pt(3)->macro_map(s, r);
      Node_pt[node_count]->x(0) = r[0];
      Node_pt[node_count]->x(1) = r[1];

      // Add the node to the boundary
      add_boundary_node(1, Node_pt[node_count]);
      // Also add to boundary 0 if it is the corner node
      if (l2 == 0)
      {
        add_boundary_node(0, Node_pt[node_count]);
      }

      // Increment the node number
      node_count++;
    }
    // Final node belongs to element 2
    // Local node number
    jnod_local = (n_p - 1) * n_p;
    unsigned jnod_local_old = (n_p - 1) * n_p;

    // Set the pointer from the element to the node
    finite_element_pt(3)->node_pt(jnod_local) =
      finite_element_pt(2)->node_pt(jnod_local_old);

    // Loop over all elements and set macro element pointer
    // to enable MacroElement-based node update
    unsigned n_element = this->nelement();
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to full element type
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(e));

      // Set pointer to macro element
      el_pt->set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }

    // Setup boundary element lookup schemes
    setup_boundary_element_info();
  }
} // namespace oomph
#endif
