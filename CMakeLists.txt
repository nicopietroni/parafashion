cmake_minimum_required(VERSION 3.1)
project(parafashion)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Qt5 COMPONENTS Widgets OpenGL REQUIRED) # redundant with libigl

# glew
include_directories(lib/glew/include)
# if this doesn't work try this https://github.com/ArthurSonzogni/OpenGL_CMake_Skeleton
#add_definitions(-DGLEW_STATIC)
#add_subdirectory(lib/glew EXCLUDE_FROM_ALL)
#target_link_libraries(main PRIVATE libglew_static)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
#    ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_definitions(-DIGL_VIEWER_VIEWER_QUIET)

add_subdirectory(lib/libigl)
add_subdirectory(lib/vcglib)

include_directories(lib/Directional/include)

# libraries which will be included in all apps:
SET(COMMON_LIBRARIES 
  igl::core
  #igl::opengl_glfw 
  #igl::png 
  #igl::opengl_glfw_imgui 
  Qt5::Widgets 
  Qt5::OpenGL 
  ${OPENGL_LIBRARIES} 
  ${GLUT_LIBRARY}
  ${GLEW_LIBRARIES}
  VGCLib
)

file(GLOB SRCFILES src/*.cpp )

add_executable(parafashion ${SRCFILES} apps/parafashion.cpp)
target_link_libraries(parafashion ${COMMON_LIBRARIES})
include_directories(include)

# TODO: other app(s), parafashion_commandline