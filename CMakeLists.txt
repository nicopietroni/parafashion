cmake_minimum_required(VERSION 3.1)
project(parafashion)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

option(IGL_STATIC_LIBRARY     "Use libIGL as static librarie" ON)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Qt5 COMPONENTS Widgets OpenGL Svg REQUIRED) # redundant with libigl ?

# glew
# see https://github.com/ArthurSonzogni/OpenGL_CMake_Skeleton
add_definitions(-DGLEW_STATIC)
add_subdirectory(lib/glew-cmake EXCLUDE_FROM_ALL)

add_subdirectory(lib/libigl)
add_subdirectory(lib/vcglib)
include_directories(lib/vcglib)
include_directories(lib/tracing/include)
include_directories(lib/AntTweakBar/include)
include_directories(lib/Clipper)

add_subdirectory(lib/json)
include_directories(lib/json/include)

include_directories(lib/libigl/include)

link_directories(AntTweakBar lib/AntTweakBar/lib)

include_directories(lib/Directional/include)

# libraries which will be included in all apps:
SET(COMMON_LIBRARIES 
  Qt5::Widgets 
  Qt5::OpenGL 
  Qt5::Svg
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARY}
  vcglib
  AntTweakBar
  libglew_static
)

file(GLOB SRCFILES src/*.cpp 
    lib/vcglib/wrap/ply/plylib.cpp # TODO cleanup ?
    lib/vcglib/wrap/gui/trackball.cpp 
    lib/vcglib/wrap/gui/trackmode.cpp
    lib/vcglib/wrap/qt/anttweakbarMapperNew.cpp
    lib/vcglib/wrap/qt/Outline2ToQImage.cpp
    lib/Clipper/clipper.cpp
)

file(GLOB SRCPARAM src/param/*.cpp)

add_executable(parafashion ${SRCFILES} apps/parafashion.cpp ${SRCPARAM})
target_link_libraries(parafashion ${COMMON_LIBRARIES} nlohmann_json) # igl::triangle) weird problems with igl::triangle...
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

file(GLOB PARAMSRC src/param/*.cpp)
add_executable(param apps/param.cpp ${PARAMSRC})
target_link_libraries(param igl::opengl_glfw igl::core igl::opengl_glfw_imgui)
target_include_directories(param PUBLIC include/param)