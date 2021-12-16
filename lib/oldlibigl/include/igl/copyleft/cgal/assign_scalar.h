// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_ASSIGN_SCALAR_H
#define IGL_COPYLEFT_CGAL_ASSIGN_SCALAR_H
#include "../../igl_inline.h"
#include "CGAL_includes.hpp"
namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Inputs:
      //   cgal  cgal scalar
      // Outputs:
      //   d  output scalar
      IGL_INLINE void assign_scalar(
        const CGAL::Epeck::FT & cgal,
        CGAL::Epeck::FT & d);
      IGL_INLINE void assign_scalar(
        const CGAL::Epeck::FT & cgal,
        double & d);
      IGL_INLINE void assign_scalar(
        const CGAL::Epeck::FT & cgal,
        float& d);
      IGL_INLINE void assign_scalar(
        const double & c,
        double & d);
      IGL_INLINE void assign_scalar(
        const float& c,
        float & d);
      IGL_INLINE void assign_scalar(
        const float& c,
        double& d);
    }
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "assign_scalar.cpp"
#endif
#endif
