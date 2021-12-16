// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "right_axis.h"
#ifndef IGL_NO_OPENGL

#include "OpenGL_convenience.h"

IGL_INLINE void igl::right_axis(double * x, double * y, double * z)
{
  double mv[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mv);
  igl::right_axis(mv,x,y,z);
}

IGL_INLINE void igl::right_axis(const double * mv,double * x, double * y, double * z)
{
  *x = -mv[0*4+0];
  *y = -mv[1*4+0];
  *z = -mv[2*4+0];
}

#endif
