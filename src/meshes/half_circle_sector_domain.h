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
// Include guards
#ifndef OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER
#define OOMPH_QUARTER_CIRCLE_SECTOR_DOMAIN_HEADER

// Generic oomph-lib includes
#include "../generic/quadtree.h"
#include "../generic/domain.h"
#include "../generic/geom_objects.h"

namespace oomph
{
  //=================================================================
  /// Circular sector as domain. Domain is bounded by
  /// curved boundary which is represented by a GeomObject. Domain is
  /// parametrised by three macro elements.
  //=================================================================
  class HalfCircleSectorDomain : public Domain
  {
  public:
    /// Constructor: Pass boundary object and start and end coordinates
    /// and fraction along boundary object where outer ring is divided.
    HalfCircleSectorDomain(GeomObject* boundary_geom_object_pt, double& x_centr)
      : Wall_pt(boundary_geom_object_pt), X_centr(x_centr)
    {
      // Calculate pi
      pi = 4.0 * atan(1.0);

      // There are three macro elements
      unsigned nmacro = 4;

      // Resize
      Macro_element_pt.resize(nmacro);

      // Create macro elements
      for (unsigned i = 0; i < nmacro; i++)
      {
        Macro_element_pt[i] = new QMacroElement<2>(this, i);
      }
    }


    /// Broken copy constructor
    HalfCircleSectorDomain(const HalfCircleSectorDomain&) = delete;

    /// Broken assignment operator
    void operator=(const HalfCircleSectorDomain&) = delete;

    /// Destructor: empty; cleanup done in base class
    ~HalfCircleSectorDomain() {}

    /// Vector representation of the  i_macro-th macro element
    /// boundary i_direct (N/S/W/E) at time level t
    /// (t=0: present; t>0: previous):
    /// f(s). Note that the local coordinate \b s is a 1D
    /// Vector rather than a scalar -- this is unavoidable because
    /// this function implements the pure virtual function in the
    /// Domain base class.
    void macro_element_boundary(const unsigned& t,
                                const unsigned& i_macro,
                                const unsigned& i_direct,
                                const Vector<double>& s,
                                Vector<double>& f);


  private:
    // pi
    double pi;

    /// Pointer to geometric object that represents the curved wall
    GeomObject* Wall_pt;

    /// x-displacement
    double X_centr;

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_bot_left_N(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_bot_left_W(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_bot_left_S(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_bot_left_E(const unsigned& t,
                      const Vector<double>& zeta,
                      Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_N(const unsigned& t,
                 const Vector<double>& zeta,
                 Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_W(const unsigned& t,
                 const Vector<double>& zeta,
                 Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_S(const unsigned& t,
                 const Vector<double>& zeta,
                 Vector<double>& f);

    /// Boundary of top left macro element zeta \f$ \in [-1,1] \f$
    void r_top_E(const unsigned& t,
                 const Vector<double>& zeta,
                 Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_N(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_W(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_S(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of bottom right macro element zeta \f$ \in [-1,1] \f$
    void r_bot_right_E(const unsigned& t,
                       const Vector<double>& zeta,
                       Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_N(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_E(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_S(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);

    /// Boundary of central box macro element zeta \f$ \in [-1,1] \f$
    void r_centr_W(const unsigned& t,
                   const Vector<double>& zeta,
                   Vector<double>& f);
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Vector representation of the  imacro-th macro element
  /// boundary idirect (N/S/W/E) at time level t (t=0: present; t>0: previous):
  /// f(s)
  //=================================================================
  void HalfCircleSectorDomain::macro_element_boundary(const unsigned& t,
                                                      const unsigned& imacro,
                                                      const unsigned& idirect,
                                                      const Vector<double>& s,
                                                      Vector<double>& f)
  {
    // Which macro element?
    // --------------------
    switch (imacro)
    {
      using namespace QuadTreeNames;

        // Macro element 0: Central box
      case 0:

        // Which direction?
        if (idirect == N)
        {
          r_centr_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_centr_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_centr_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_centr_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 1: Bottom right
      case 1:

        // Which direction?
        if (idirect == N)
        {
          r_bot_right_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_bot_right_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_bot_right_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_bot_right_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 2:Top
      case 2:

        // Which direction?
        if (idirect == N)
        {
          r_top_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_top_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_top_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_top_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

        // Macro element 3:Top
      case 3:

        // Which direction?
        if (idirect == N)
        {
          r_bot_left_N(t, s, f);
        }
        else if (idirect == S)
        {
          r_bot_left_S(t, s, f);
        }
        else if (idirect == W)
        {
          r_bot_left_W(t, s, f);
        }
        else if (idirect == E)
        {
          r_bot_left_E(t, s, f);
        }
        else
        {
          std::ostringstream error_stream;
          error_stream << "idirect is " << idirect << " not one of N, S, E, W"
                       << std::endl;

          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        break;

      default:

        // Error
        std::ostringstream error_stream;
        error_stream << "Wrong imacro " << imacro << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //=================================================================
  /// Northern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void HalfCircleSectorDomain::r_top_N(const unsigned& t,
                                       const Vector<double>& s,
                                       Vector<double>& f)
  {
    Vector<double> x(1);

    // Coordinate along wall
    x[0] = 3.0 * pi / 4.0 + (pi / 4.0 - 3.0 * pi / 4.0) * 0.5 * (s[0] + 1.0);

    Wall_pt->position(t, x, f);

    f[0] += X_centr;
  }


  //=================================================================
  /// Western edge of top left macro element \f$s \in [-1,1] \f$
  //=================================================================
  void HalfCircleSectorDomain::r_top_W(const unsigned& t,
                                       const Vector<double>& s,
                                       Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_top(2);
    x[0] = 3.0 * pi / 4.0;
    Wall_pt->position(t, x, r_top);

    Vector<double> r_bottom(2, 0.0);
    x[0] = pi;
    Wall_pt->position(t, x, r_bottom);
    r_bottom[0] = r_bottom[0] / 2.0;

    x[0] = pi / 2.0;
    Vector<double> r_tmp(2, 0.0);
    Wall_pt->position(t, x, r_tmp);
    r_bottom[1] = r_tmp[1] / 2.0;

    for (unsigned j = 0; j < 2; j++)
    {
      f[j] = r_bottom[j] + 0.5 * (s[0] + 1.0) * (r_top[j] - r_bottom[j]);
    }

    f[0] += X_centr;
  }


  //=================================================================
  /// Southern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void HalfCircleSectorDomain::r_top_S(const unsigned& t,
                                       const Vector<double>& s,
                                       Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_top(2);
    x[0] = pi / 2.0;
    Wall_pt->position(t, x, r_top);
    r_top[1] = r_top[1] / 2.0;

    Vector<double> r_left(2);
    x[0] = pi;
    Wall_pt->position(t, x, r_left);
    r_left[0] = r_left[0] / 2.0;

    f[0] = r_left[0] + 0.5 * (1.0 + s[0]) * (-r_left[0] - r_left[0]);
    f[1] = r_top[1];

    f[0] += X_centr;
  }


  //=================================================================
  /// Eastern edge of top left macro element \f$ s \in [-1,1] \f$
  //=================================================================
  void HalfCircleSectorDomain::r_top_E(const unsigned& t,
                                       const Vector<double>& s,
                                       Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_top(2);
    x[0] = pi / 4.0;
    Wall_pt->position(t, x, r_top);

    Vector<double> r_bottom(2, 0.0);
    x[0] = 0.0;
    Wall_pt->position(t, x, r_bottom);
    r_bottom[0] = r_bottom[0] / 2.0;

    x[0] = pi / 2.0;
    Vector<double> r_tmp(2, 0.0);
    Wall_pt->position(t, x, r_tmp);
    r_bottom[1] = r_tmp[1] / 2.0;

    for (unsigned j = 0; j < 2; j++)
    {
      f[j] = r_bottom[j] + 0.5 * (s[0] + 1.0) * (r_top[j] - r_bottom[j]);
    }

    f[0] += X_centr;
  }


  //=================================================================
  /// Northern edge of bottom right macro element
  //=================================================================
  void HalfCircleSectorDomain::r_bot_right_N(const unsigned& t,
                                             const Vector<double>& s,
                                             Vector<double>& f)
  {
    r_top_E(t, s, f);
  }

  //=================================================================
  /// Western edge of bottom right macro element
  //=================================================================
  void HalfCircleSectorDomain::r_bot_right_W(const unsigned& t,
                                             const Vector<double>& s,
                                             Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_right(2);
    x[0] = 0.0;
    Wall_pt->position(t, x, r_right);
    r_right[0] = r_right[0] / 2.0;

    Vector<double> r_top(2, 0.0);
    x[0] = pi / 2.0;
    Wall_pt->position(t, x, r_top);
    r_top[1] = r_top[1] / 2.0;

    f[0] = r_right[0];
    f[1] = 0.5 * (1.0 + s[0]) * r_top[1];

    f[0] += X_centr;
  }

  //=================================================================
  /// Southern edge of bottom right macro element
  //=================================================================
  void HalfCircleSectorDomain::r_bot_right_S(const unsigned& t,
                                             const Vector<double>& s,
                                             Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_left(2);
    x[0] = 0.0;
    Wall_pt->position(t, x, r_left);
    r_left[0] = r_left[0] / 2.0;

    f[0] = r_left[0] + 0.5 * (1.0 + s[0]) * (2.0 * r_left[0] - r_left[0]);
    f[1] = 0.0;

    f[0] += X_centr;
  }

  //=================================================================
  /// Eastern edge of bottom right macro element
  //=================================================================
  void HalfCircleSectorDomain::r_bot_right_E(const unsigned& t,
                                             const Vector<double>& s,
                                             Vector<double>& f)
  {
    Vector<double> x(1);

    // Coordinate along wall
    x[0] = 0.5 * (1.0 + s[0]) * pi / 4.0;

    Wall_pt->position(t, x, f);

    f[0] += X_centr;
  }


  //=================================================================
  /// Northern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_centr_N(const unsigned& t,
                                         const Vector<double>& s,
                                         Vector<double>& f)
  {
    r_top_S(t, s, f);
  }


  //=================================================================
  /// Eastern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_centr_E(const unsigned& t,
                                         const Vector<double>& s,
                                         Vector<double>& f)
  {
    r_bot_right_W(t, s, f);
  }


  //=================================================================
  /// Southern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_centr_S(const unsigned& t,
                                         const Vector<double>& s,
                                         Vector<double>& f)
  {
    Vector<double> x(1);

    // Bottom right corner
    Vector<double> r_left(2);
    x[0] = pi;
    Wall_pt->position(t, x, r_left);
    r_left[0] = r_left[0] / 2.0;

    f[0] = r_left[0] + 0.5 * (s[0] + 1.0) * (-r_left[0] - r_left[0]);
    f[1] = 0.0;

    f[0] += X_centr;
  }


  //=================================================================
  /// Western  edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_centr_W(const unsigned& t,
                                         const Vector<double>& s,
                                         Vector<double>& f)
  {
    Vector<double> x(1);

    Vector<double> r_left(2);
    x[0] = pi;
    Wall_pt->position(t, x, r_left);
    r_left[0] = r_left[0] / 2.0;

    Vector<double> r_top(2);
    x[0] = pi / 2.0;
    Wall_pt->position(t, x, r_top);
    r_top[1] = r_top[1] / 2.0;

    f[0] = r_left[0];
    f[1] = 0.5 * (1.0 + s[0]) * r_top[1];

    f[0] += X_centr;
  }

  //=================================================================
  /// Northern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_bot_left_N(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    Vector<double> s_tmp(1, -s[0]);
    r_top_W(t, s, f);
  }


  //=================================================================
  /// Eastern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_bot_left_E(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    r_centr_W(t, s, f);
  }


  //=================================================================
  /// Southern edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_bot_left_S(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    Vector<double> x(1);

    // Bottom right corner
    Vector<double> r_left(2);
    x[0] = pi;
    Wall_pt->position(t, x, r_left);

    f[0] = 2.0 * r_left[0] + 0.5 * (s[0] + 1.0) * (r_left[0] - 2.0 * r_left[0]);
    f[1] = 0.0;

    f[0] += X_centr;
  }


  //=================================================================
  /// Western  edge of central box
  //=================================================================
  void HalfCircleSectorDomain::r_bot_left_W(const unsigned& t,
                                            const Vector<double>& s,
                                            Vector<double>& f)
  {
    Vector<double> x(1);
    x[0] = 3.0 * pi / 4.0 + 0.5 * (1.0 - s[0]) * (pi - 3.0 * pi / 4.0);
    Wall_pt->position(t, x, f);

    f[0] += X_centr;
  }


} // namespace oomph

#endif
