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

if (NOT TARGET igl::core)
  add_subdirectory(lib/libigl)
endif() 
add_subdirectory(lib/vcglib)
add_subdirectory(lib/garment-flattening)
add_subdirectory(lib/AntTweakBar)
include_directories(lib/vcglib)
include_directories(lib/tracing/include)
include_directories(lib/AntTweakBar/TweakBar/include)
include_directories(lib/Clipper)

if (NOT TARGET nlohmann_json)
  add_subdirectory(lib/json)
endif() 
include_directories(lib/json/include)
include_directories(lib/libigl/include)

link_directories(AntTweakBar lib/AntTweakBar/lib)

include_directories(lib/Directional/include)

# libraries which will be included in all apps:
SET(VCG_UI_LIBS 
  Qt5::Widgets 
  Qt5::OpenGL 
  Qt5::Svg
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARY}
  vcglib
  AntTweakBar
  libglew_static
)

file(GLOB SRC_VCGLIB src/*.cpp 
    lib/vcglib/wrap/ply/plylib.cpp
    lib/vcglib/wrap/gui/trackball.cpp 
    lib/vcglib/wrap/gui/trackmode.cpp
    lib/vcglib/wrap/qt/anttweakbarMapperNew.cpp
    lib/vcglib/wrap/qt/Outline2ToQImage.cpp
    lib/Clipper/clipper.cpp
)

file(GLOB SRCPARAM src/param/*.cpp)

add_executable(parafashion ${SRC_VCGLIB} apps/parafashion.cpp ${SRCPARAM})
target_link_libraries(parafashion ${VCG_UI_LIBS} nlohmann_json woven_param)
include_directories(include)

add_executable(param apps/param.cpp ${SRCPARAM})
target_link_libraries(param igl::opengl_glfw igl::core igl::opengl_glfw_imgui)
target_include_directories(param PUBLIC include/param)
