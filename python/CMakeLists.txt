find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(sayram2d_py sayram2d_py.cc)
target_link_libraries(sayram2d_py PRIVATE sayram2d)
set_target_properties(sayram2d_py PROPERTIES
  OUTPUT_NAME sayram2d
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
