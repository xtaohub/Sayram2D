cmake_minimum_required(VERSION 3.16)
project(sayram2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

find_package(Eigen3 3.4 QUIET CONFIG)
if(NOT TARGET Eigen3::Eigen)
  find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Dense
    HINTS ENV EIGEN3_ROOT
      /opt/homebrew/include/eigen3
      /opt/homebrew/opt/eigen/include/eigen3
      /usr/local/include/eigen3
      /usr/include/eigen3)
  if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR
      "Eigen3 not found. Set Eigen3_DIR or EIGEN3_INCLUDE_DIR.")
  endif()
endif()

add_library(sayram2d
  src/mesh.cc
  src/boundary_conditions.cc
  src/solver.cc
)
target_include_directories(sayram2d PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src/include
)
if(TARGET Eigen3::Eigen)
  target_link_libraries(sayram2d PUBLIC Eigen3::Eigen)
else()
  target_include_directories(sayram2d SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()
target_compile_definitions(sayram2d PUBLIC EIGEN_NO_DEBUG)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(sayram2d PRIVATE -Wall -Wextra)
endif()

add_executable(sayram2d_example1
  src/examples/example1.cc
)
target_link_libraries(sayram2d_example1 PRIVATE sayram2d)

add_executable(sayram2d_example2
  src/examples/example2.cc
)
target_link_libraries(sayram2d_example2 PRIVATE sayram2d)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(sayram2d_example1 PRIVATE -Wall -Wextra)
  target_compile_options(sayram2d_example2 PRIVATE -Wall -Wextra)
endif()

option(SAYRAM2D_BUILD_APP "Build the standalone HDF5 application" ON)
if(SAYRAM2D_BUILD_APP)
  # HDF5 is an application-only dependency. Configure with
  # -DSAYRAM2D_BUILD_APP=OFF for a library-only build without HDF5.
  find_package(HDF5 COMPONENTS C QUIET)
  if(NOT HDF5_FOUND)
    find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
      HINTS ENV HDF5_ROOT /opt/homebrew/include /opt/homebrew/opt/hdf5/include
      /usr/local/include)
    find_library(HDF5_LIBRARY NAMES hdf5
      HINTS ENV HDF5_ROOT /opt/homebrew/lib /opt/homebrew/opt/hdf5/lib
      /usr/local/lib)
    if(NOT HDF5_INCLUDE_DIR OR NOT HDF5_LIBRARY)
      message(FATAL_ERROR
        "HDF5 not found. Set HDF5_ROOT or configure with "
        "-DSAYRAM2D_BUILD_APP=OFF for a library-only build.")
    endif()
  endif()

  add_library(sayram2d_app_core STATIC
    src/app/config.cc
    src/app/io.cc
    src/examples/albert_young/albert_young.cc
  )
  target_include_directories(sayram2d_app_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/app
    ${CMAKE_CURRENT_SOURCE_DIR}/src/examples/albert_young
  )
  target_link_libraries(sayram2d_app_core PUBLIC sayram2d)
  if(TARGET HDF5::HDF5)
    target_link_libraries(sayram2d_app_core PUBLIC HDF5::HDF5)
  else()
    target_include_directories(sayram2d_app_core PUBLIC ${HDF5_INCLUDE_DIR})
    target_link_libraries(sayram2d_app_core PUBLIC ${HDF5_LIBRARY})
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sayram2d_app_core PRIVATE -Wall -Wextra)
  endif()

  add_executable(sayram2d_cli src/examples/albert_young/main.cc)
  target_link_libraries(sayram2d_cli PRIVATE sayram2d_app_core)
  set_target_properties(sayram2d_cli PROPERTIES
    OUTPUT_NAME sayram2d
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sayram2d_cli PRIVATE -Wall -Wextra)
  endif()
endif()

option(SAYRAM2D_BUILD_PYTHON "Build the Python bindings" OFF)
if(SAYRAM2D_BUILD_PYTHON)
  add_subdirectory(python)
endif()

if(BUILD_TESTING AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sayram2d.cc)
  add_executable(sayram2d_tests tests/test_sayram2d.cc)
  target_link_libraries(sayram2d_tests PRIVATE sayram2d)
  add_test(NAME sayram2d_unit COMMAND sayram2d_tests)
  set_tests_properties(sayram2d_unit PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(sayram2d_reference tests/new_reference.cc)
  target_link_libraries(sayram2d_reference PRIVATE sayram2d)
  add_test(NAME sayram2d_reference_generate
    COMMAND sayram2d_reference build/new_reference)
  set_tests_properties(sayram2d_reference_generate PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(sayram2d_albert_young tests/albert_young_reference.cc)
  target_link_libraries(sayram2d_albert_young PRIVATE sayram2d)
  add_test(NAME sayram2d_albert_young_run
    COMMAND sayram2d_albert_young ${CMAKE_CURRENT_SOURCE_DIR}/data/albert_young/ALBERT_YOUNG_80x80.bin
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/albert_young/new_reference.csv)
  set_tests_properties(sayram2d_albert_young_run PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME sayram2d_reference_compare
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_legacy.py
        ${CMAKE_CURRENT_SOURCE_DIR}/build/new_reference
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/legacy/reference)
    set_tests_properties(sayram2d_reference_compare PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS sayram2d_reference_generate)

    add_test(NAME sayram2d_albert_young_compare
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_albert_young.py
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/albert_young/legacy_reference.csv
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/albert_young/new_reference.csv)
    set_tests_properties(sayram2d_albert_young_compare PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS sayram2d_albert_young_run)
  endif()

  if(TARGET sayram2d_app_core)
    add_executable(sayram2d_app_tests tests/test_app.cc)
    target_link_libraries(sayram2d_app_tests PRIVATE sayram2d_app_core)
    add_test(NAME sayram2d_app_unit COMMAND sayram2d_app_tests)
    set_tests_properties(sayram2d_app_unit PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME sayram2d_app_run
      COMMAND $<TARGET_FILE:sayram2d_cli>
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/app_test.ini)
    set_tests_properties(sayram2d_app_run PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    if(Python3_Interpreter_FOUND)
      add_test(NAME sayram2d_app_compare
        COMMAND ${Python3_EXECUTABLE}
          ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_app_hdf5.py
          ${CMAKE_CURRENT_SOURCE_DIR}/build/app_test_output/AppTest/AppTest_data.h5
          ${CMAKE_CURRENT_SOURCE_DIR}/tests/albert_young/legacy_reference.csv)
      set_tests_properties(sayram2d_app_compare PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS sayram2d_app_run)
    endif()
  endif()
endif()
