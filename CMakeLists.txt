cmake_minimum_required(VERSION 3.16)
project(sayram2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Dense REQUIRED)

add_library(sayram2d
  src/mesh.cpp
  src/boundary_conditions.cpp
  src/solver.cpp
)
target_include_directories(sayram2d PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_include_directories(sayram2d SYSTEM PUBLIC
  ${EIGEN3_INCLUDE_DIR}
)
target_compile_definitions(sayram2d PUBLIC EIGEN_NO_DEBUG)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(sayram2d PRIVATE -Wall -Wextra)
endif()

add_executable(sayram2d_example
  src/main.cc
)
target_link_libraries(sayram2d_example PRIVATE sayram2d)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(sayram2d_example PRIVATE -Wall -Wextra)
endif()
