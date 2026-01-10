cmake_minimum_required(VERSION 3.12) # Bumped slightly for better GLOB handling
project(ROS_Solver VERSION 1.0 LANGUAGES Fortran)

# 1. Prevent in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In-source builds are not allowed. Please use: cmake -B build")
endif()

# 2. Module & Compiler Setup
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Add flags (removed -fcheck=all for Release builds if desired, but kept as requested)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -g -Wall -fcheck=all -fbacktrace")

# 3. Dependencies
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

# 4. Core Object Library
# ---------------------------------------------------------------------------
# Grab ALL .f90 files in src recursively
file(GLOB_RECURSE ALL_SOURCES "src/*.f90")

# Create a clean list for the library modules
set(LIB_SOURCES ${ALL_SOURCES})

# EXCLUDE src/main.f90 (It is the executable, not part of the library)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*/src/main\\.f90$")

# EXCLUDE src/tests/* (They are standalone test programs)
list(FILTER LIB_SOURCES EXCLUDE REGEX ".*/src/tests/.*")

# Create the library with the auto-detected list
add_library(ros_core OBJECT ${LIB_SOURCES})
# ---------------------------------------------------------------------------

# 5. Main Executable
add_executable(ros_main src/main.f90)
target_sources(ros_main PUBLIC $<TARGET_OBJECTS:ros_core>)
target_link_libraries(ros_main PUBLIC LAPACK::LAPACK OpenMP::OpenMP_Fortran)

# 6. Test Suite
enable_testing()

# Define the macro for creating tests
macro(create_test test_name source_file)
    add_executable(${test_name} ${source_file})
    target_sources(${test_name} PUBLIC $<TARGET_OBJECTS:ros_core>)
    target_link_libraries(${test_name} PUBLIC LAPACK::LAPACK OpenMP::OpenMP_Fortran)
    set_target_properties(${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
    add_test(NAME ${test_name} COMMAND ${test_name})
endmacro()

# ---------------------------------------------------------------------------
# Automatically discover and create tests from src/tests/*.f90
file(GLOB TEST_SOURCES "src/tests/*.f90")

foreach(test_src ${TEST_SOURCES})
    # Extract filename without extension (e.g., "src/tests/test_atmosphere.f90" -> "test_atmosphere")
    get_filename_component(test_name ${test_src} NAME_WE)
    
    # Generate the test target
    create_test(${test_name} ${test_src})
endforeach()
# ---------------------------------------------------------------------------