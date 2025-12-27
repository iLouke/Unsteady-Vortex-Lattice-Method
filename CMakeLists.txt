cmake_minimum_required(VERSION 3.10)
project(MyVLM VERSION 1.0 LANGUAGES Fortran)

# -----------------------------------------------------------------------------
# 1. Pre-build Checks
# -----------------------------------------------------------------------------
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR "In-source builds are not allowed. Please run cmake from a 'build' folder.")
endif()

# Store .mod files in a separate directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# -----------------------------------------------------------------------------
# 2. Compiler Flags & Dependencies
# -----------------------------------------------------------------------------

# Your specific compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -g -Wall -fcheck=all -fbacktrace")

# --- FIND LAPACK ---
find_package(LAPACK REQUIRED)
message(STATUS "LAPACK Found: ${LAPACK_FOUND}")

# --- FIND OPENMP ---
find_package(OpenMP REQUIRED)

if(OpenMP_Fortran_FOUND)
    message(STATUS "OpenMP Found. Multithreading enabled.")
    # It is cleaner in modern CMake to link the target later, 
    # but adding flags here ensures compilation uses them during object creation.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(WARNING "OpenMP not found. Running in serial mode.")
endif()

# -----------------------------------------------------------------------------
# 3. Source Definitions (Updated for new folder structure)
# -----------------------------------------------------------------------------

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
    src/solvers/mod_integration.f90
    src/solvers/mod_linalg.f90
    src/solvers/mod_vector.f90
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

set(MAIN_SOURCE
    src/main.f90
)

# -----------------------------------------------------------------------------
# 4. Build Targets
# -----------------------------------------------------------------------------

add_executable(vlm_solver
    ${CORE_SOURCES}
    ${GEO_SOURCES}
    ${MAIN_SOURCE}
)

# -----------------------------------------------------------------------------
# 5. Link Libraries
# -----------------------------------------------------------------------------

# Link LAPACK
target_link_libraries(vlm_solver ${LAPACK_LIBRARIES})

# Link OpenMP
if(OpenMP_Fortran_FOUND)
    target_link_libraries(vlm_solver OpenMP::OpenMP_Fortran)
endif()