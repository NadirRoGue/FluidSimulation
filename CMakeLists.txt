cmake_minimum_required(VERSION 2.8)

project(practica3)
add_executable(practica3
  src/index2.cpp
  src/fluid2Exercise.cpp
  src/fluidVisualizer2.cpp
  src/main.cpp
  src/scene.cpp
  src/vec2.cpp)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pcgsolver)

target_link_libraries(practica3 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
