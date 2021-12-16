// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_UNIFORM_TYPE_TO_STRING_H
#define IGL_UNIFORM_TYPE_TO_STRING_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include <string>

#include "OpenGL_convenience.h"

namespace igl
{
  // Convert a GL uniform variable type (say, returned from
  // glGetActiveUniform) and output a string naming that type
  // Inputs:
  //   type  enum for given type
  // Returns string name of that type
  IGL_INLINE std::string uniform_type_to_string(const GLenum type);
}

#ifndef IGL_STATIC_LIBRARY
#  include "uniform_type_to_string.cpp"
#endif

#endif
#endif
