cmake_minimum_required(VERSION 3.10)
project(MyVLM VERSION 1.0 LANGUAGES Fortran)

# -----------------------------------------------------------------------------
# 1. Pre-build Checks & Settings
# -----------------------------------------------------------------------------
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR "In-source builds not allowed. Run cmake from a 'build' folder.")
endif()

# Enable Testing (Standard CTest functionality)
enable_testing()

# Module Directory setup
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
# IMPORTANT: Ensure all targets can find the .mod files generated here
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Compiler Flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -g -Wall -fcheck=all -fbacktrace")

# --- FIND DEPENDENCIES ---
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

if(OpenMP_Fortran_FOUND)
    message(STATUS "OpenMP Found. Multithreading enabled.")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(WARNING "OpenMP not found. Running in serial mode.")
endif()

# -----------------------------------------------------------------------------
# 2. Source Definitions
# -----------------------------------------------------------------------------
# (Your existing lists remain exactly the same)

set(CORE_SOURCES
    src/core/mod_constants.f90
    src/core/mod_precision.f90
)

set(GEO_SOURCES
    src/geometry/mod_geom_tools.f90
    src/geometry/mod_mesh_types.f90
)

set(IO_SOURCES
    src/io/mod_nastran.f90
    src/io/mod_vtk.f90
)

set(MATH_SOURCES
    src/math/mod_integration.f90
    src/math/mod_linalg.f90
    src/math/mod_analytic_geometry.f90
)

set(PHYSICS_SOURCES
    src/physics/mod_atmosphere.f90
)

set(SINGULARITIES_SOURCES
    src/singularities/mod_doublet.f90
    src/singularities/mod_potential.f90
    src/singularities/mod_source.f90
    src/singularities/mod_vortex.f90
)

set(SOLVERS_SOURCES
    src/solvers/mod_panel_solver.f90
    src/solvers/mod_vlm_solver.f90
)

# Note: Main is NOT in the library lists
set(MAIN_SOURCE src/main.f90)

# -----------------------------------------------------------------------------
# 3. Create the Core Library
# -----------------------------------------------------------------------------
# This compiles all your modules ONCE into a static library.
add_library(vlm_core STATIC
    ${CORE_SOURCES}
    # ${GEO_SOURCES}
    # ${IO_SOURCES}
    ${MATH_SOURCES}
    ${PHYSICS_SOURCES}
    # ${SINGULARITIES_SOURCES}
    # ${SOLVERS_SOURCES}
)

# Link external dependencies to the core library so they propagate
target_link_libraries(vlm_core PUBLIC ${LAPACK_LIBRARIES})
if(OpenMP_Fortran_FOUND)
    target_link_libraries(vlm_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

# -----------------------------------------------------------------------------
# 4. Main Executable
# -----------------------------------------------------------------------------
add_executable(vlm_solver ${MAIN_SOURCE})
# Link to your own library
target_link_libraries(vlm_solver PRIVATE vlm_core)

# -----------------------------------------------------------------------------
# 5. Tests
# -----------------------------------------------------------------------------
message(STATUS "Configuring Tests...")

# --- Test 1: Atmosphere Model ---
add_executable(test_atmosphere src/tests/test_atmosphere.f90)

target_link_libraries(test_atmosphere PRIVATE vlm_core)

# Register the test with CTest
add_test(NAME AtmosphereTest COMMAND test_atmosphere)