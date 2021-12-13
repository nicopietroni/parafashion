cmake_minimum_required(VERSION 3.1)
project(parafashion)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

option(IGL_STATIC_LIBRARY "Use libIGL as static library" OFF)

option(VCG_HEADER_ONLY "Use VCG as static library" ON)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Qt5 COMPONENTS Widgets OpenGL Svg Xml REQUIRED)

# glew
# see https://github.com/ArthurSonzogni/OpenGL_CMake_Skeleton
add_definitions(-DGLEW_STATIC)
add_subdirectory(lib/glew-cmake EXCLUDE_FROM_ALL)

#add_subdirectory(lib/libigl) # TODO use this instead of hard copy in field_tracing
#add_subdirectory(lib/Directional/tutorials) # TODO use this instead of hard copy in field_tracing
#add_subdirectory(lib/vcglib)
include_directories(lib/vcglib)
include_directories(lib/AntTweakBar/include)
include_directories(lib/Clipper)

link_directories(AntTweakBar lib/AntTweakBar/lib)

add_library(igl_wrap SHARED 
  lib/field_tracing/lib/libigl/include/igl/principal_curvature.h 
  lib/field_tracing/lib/libigl/include/igl/copyleft/comiso/nrosy.h)
set_target_properties(igl_wrap PROPERTIES LINKER_LANGUAGE CXX)

# ------------

include_directories(lib/field_tracing/lib/libigl/external/nanogui/ext/eigen) # TODO replace this with libigl's Eigen
include_directories(lib/field_tracing) # TODO update
include_directories(lib/field_tracing/tracing)
include_directories(lib/field_tracing/lib/libigl/include) # TODO REMOVE


#add_subdirectory(lib/field_tracing/lib/libigl/external/nanogui/ext/eigen/Eigen)
#add_subdirectory(lib/field_tracing)

#add_library(field_tracing SHARED lib/field_tracing/main.cpp lib/field_tracing/main.cpp)
#target_include_directories(field_tracing PUBLIC
#  lib/field_tracing
#  lib/field_tracing/tracing
#  lib/field_tracing/lib/libigl/external/nanogui/ext/eigen
#  lib/field_tracing/lib/libigl/include
#)
#target_link_libraries(field_tracing
#  vcg_wrap
#  Qt5::Widgets 
#  Qt5::OpenGL 
#  Qt5::Svg
#  libglew_static
#)
# ------------

#add_library(vcg_wrap SHARED 
#  lib/vcglib/wrap/ply/plylib.cpp # TODO cleanup ?
#  lib/vcglib/wrap/gui/trackball.cpp 
#  lib/vcglib/wrap/gui/trackmode.cpp
#  lib/vcglib/wrap/qt/anttweakbarMapperNew.cpp
#  lib/vcglib/wrap/qt/Outline2ToQImage.cpp
#s)
#target_include_directories(vcg_wrap PUBLIC
#  lib/vcglib/
#)
#target_link_libraries(vcg_wrapfield_tracing
#  igl::core
#  igl_wrap
#  #tutorials
#  Qt5::Widgets 
#  Qt5::OpenGL 
#  Qt5::Svg
#  libglew_static
#)

#include_directories(lib/Directional/include)



#add_library(vcg_wrap SHARED ${SOURCES_WRAP} ${VCG_HEADERS})
#target_include_directories(vcg_wrap PUBLIC  lib/vcg/ # TODO /include directory lib/AntTweakBar/include	lib/vcg/eigenlib/) 
#target_link_libraries(vcg_wrap PUBLIC	Qt5::OpenGL 	Qt5::Widgets 	Qt5::Svg	Qt5::Xml	AntTweakBar	${OPENGL_LIBRARIES} 	${GLUT_LIBRARY}	libglew_static)

#file(GLOB_RECURSE VCG_HEADERS "lib/vcglib/*.h")
#file(GLOB_RECURSE VCG_SOURCES "lib/vcglib/*.cpp")
#add_library(vcg_wrap ${VCG_SOURCES} ${VCG_HEADERS})
#include_directories(vcg_wrap lib/vcglib/)
#target_link_libraries(vcg_wrap PUBLIC Qt5 Qt5::Widgets Qt5::OpenGL Qt5::Svg ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} AntTweakBar libglew_static)



# libraries which will be included in all apps:
SET(COMMON_LIBRARIES 
  igl_wrap
  #field_tracing
  #vcg_wrap
  Qt5::Widgets 
  Qt5::OpenGL 
  Qt5::Svg
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARY}
  AntTweakBar
  libglew_static
)

file(GLOB SRCFILES src/*.cpp 
  lib/Clipper/clipper.cpp
  lib/vcglib/wrap/ply/plylib.cpp # TODO cleanup ?
  lib/vcglib/wrap/gui/trackball.cpp 
  lib/vcglib/wrap/gui/trackmode.cpp
  lib/vcglib/wrap/qt/anttweakbarMapperNew.cpp
  lib/vcglib/wrap/qt/Outline2ToQImage.cpp
)

add_executable(parafashion ${SRCFILES} apps/parafashion_app.cpp)
target_link_libraries(parafashion ${COMMON_LIBRARIES})
include_directories(include)

file(GLOB TEXTILESRC # TODO REMOVE...
  src/textile_param/ourglwidget.cpp 
  src/textile_param/*.cpp 
  lib/vcglib/wrap/ply/plylib.cpp
  lib/vcglib/wrap/gui/trackball.cpp
  lib/vcglib/wrap/gui/trackmode.cpp
  lib/vcglib/wrap/qt/anttweakbarMapperNew.cpp
  #lib/libigl/principal_curvature.cpp
)

add_executable(textile_param apps/textile_param.cpp ${TEXTILESRC})
target_link_libraries(textile_param 
  Qt5::Widgets 
  Qt5::OpenGL 
  Qt5::Svg
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARY}
  AntTweakBar
  libglew_static
  ${COMMON_LIBRARIES}
)
target_include_directories(textile_param PUBLIC include/textile_param)

# TODO: other app(s), parafashion_commandline