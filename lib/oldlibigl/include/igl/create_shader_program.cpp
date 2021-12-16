// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "create_shader_program.h"
#ifndef IGL_NO_OPENGL

#include "load_shader.h"
#include "print_program_info_log.h"
#include <iostream>
#include <cstdio>

IGL_INLINE bool igl::create_shader_program(
  const std::string & geom_source,
  const std::string & vert_source,
  const std::string & frag_source,
  const std::map<std::string,GLuint> & attrib,
  GLuint & id)
{
  using namespace std;
  if(vert_source == "" && frag_source == "")
  {
    cerr<<
      "create_shader_program() could not create shader program,"
      " both .vert and .frag source given were empty"<<endl;
    return false;
  }

  // create program
  id = glCreateProgram();
  if(id == 0)
  {
    cerr<<"create_shader_program() could not create shader program."<<endl;
    return false;
  }

  if(geom_source != "")
  {
    // load vertex shader
    GLuint g = igl::load_shader(geom_source.c_str(),GL_GEOMETRY_SHADER_EXT);
    if(g == 0)
    {
      cerr<<"geometry shader failed to compile."<<endl;
      return false;
    }
    glAttachShader(id,g);
  }

  if(vert_source != "")
  {
    // load vertex shader
    GLuint v = igl::load_shader(vert_source.c_str(),GL_VERTEX_SHADER);
    if(v == 0)
    {
      cerr<<"vertex shader failed to compile."<<endl;
      return false;
    }
    glAttachShader(id,v);
  }

  if(frag_source != "")
  {
    // load fragment shader
    GLuint f = igl::load_shader(frag_source.c_str(),GL_FRAGMENT_SHADER);
    if(f == 0)
    {
      cerr<<"fragment shader failed to compile."<<endl;
      return false;
    }
    glAttachShader(id,f);
  }

  // loop over attributes
  for(
    std::map<std::string,GLuint>::const_iterator ait = attrib.begin();
    ait != attrib.end();
    ait++)
  {
    glBindAttribLocation(
      id,
      (*ait).second,
      (*ait).first.c_str());
  }
  // Link program
  glLinkProgram(id);

  // print log if any
  igl::print_program_info_log(id);

  return true;
}

IGL_INLINE bool igl::create_shader_program(
  const std::string & vert_source,
  const std::string & frag_source,
  const std::map<std::string,GLuint> & attrib,
  GLuint & prog_id)
{
  return create_shader_program("",vert_source,frag_source,attrib,prog_id);
}


IGL_INLINE GLuint igl::create_shader_program(
  const std::string & geom_source,
  const std::string & vert_source,
  const std::string & frag_source,
  const std::map<std::string,GLuint> & attrib)
{
  GLuint prog_id = 0;
  create_shader_program(geom_source,vert_source,frag_source,attrib,prog_id);
  return prog_id;
}

IGL_INLINE GLuint igl::create_shader_program(
  const std::string & vert_source,
  const std::string & frag_source,
  const std::map<std::string,GLuint> & attrib)
{
  GLuint prog_id = 0;
  create_shader_program(vert_source,frag_source,attrib,prog_id);
  return prog_id;
}

#endif
