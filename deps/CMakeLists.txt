include(ExternalProject)

add_subdirectory(${PROJECT_SOURCE_DIR}/deps/fmtlib EXCLUDE_FROM_ALL)
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/pugixml)

if(EXISTS ${PROJECT_SOURCE_DIR}/deps/pybind11/CMakeLists.txt)
    add_subdirectory(${PROJECT_SOURCE_DIR}/deps/pybind11)
else()
    find_package(pybind11 CONFIG REQUIRED)
endif()


if(NOT USE_SYSTEM_EIGEN3)
ExternalProject_Add(eigen3

  SOURCE_DIR "${PROJECT_SOURCE_DIR}/deps/eigen3"
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
)
endif()

if(BUILD_QC)
endif(BUILD_QC)
