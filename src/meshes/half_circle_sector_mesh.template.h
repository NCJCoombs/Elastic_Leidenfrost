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
#ifndef OOMPH_HALF_CIRCLE_SECTOR_MESH_HEADER
#define OOMPH_HALF_CIRCLE_SECTOR_MESH_HEADER

#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/quad_mesh.h"
#include "../generic/macro_element_node_update_element.h"

// Include the headers file for domain
#include "half_circle_sector_domain.h"


namespace oomph
{
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  class GeomObject;

  template<class ELEMENT>
  class HalfCircleSectorMesh : public virtual QuadMeshBase
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh)
    HalfCircleSectorMesh(
      GeomObject* wall_pt,
      double& x_centr,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor:
    virtual ~HalfCircleSectorMesh() {}

    /// Access function to GeomObject representing wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }

    /// Access function to domain
    HalfCircleSectorDomain* domain_pt()
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to Domain
    HalfCircleSectorDomain* Domain_pt;

    /// Pointer to the geometric object that represents the curved wall
    /// (mesh boundary 1)
    GeomObject* Wall_pt;

    /// x-displacement
    double X_centr;
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////

  template<class ELEMENT>
  class RefineableHalfCircleSectorMesh
    : public HalfCircleSectorMesh<ELEMENT>,
      public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh).
    /// Adds refinement data to elements of HalfCircleSectorMesh.
    RefineableHalfCircleSectorMesh(
      GeomObject* wall_pt,
      double& x_centr,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : HalfCircleSectorMesh<ELEMENT>(wall_pt, x_centr, time_stepper_pt)
    {
      // Basic mesh has been built -- just need to setup the
      // adaptivity information:

      // Setup quadtree forest
      this->setup_quadtree_forest();
    }

    /// Destructor: Empty
    virtual ~RefineableHalfCircleSectorMesh() {}
  };
} // namespace oomph

#endif
