cmake_minimum_required(VERSION 3.12)
project(AGPL_Fortran LANGUAGES Fortran)

# Export compile commands for tools like language servers, IDEs, etc.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set Fortran compiler flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_Fortran_FLAGS "-std=f2008 -fPIC")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -Wall")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")

# Find OpenMP (used by integration module)
find_package(OpenMP REQUIRED)

# Set module output directory
set(FORTRAN_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Include directories for modules
include_directories(${FORTRAN_MODULE_DIRECTORY})

# Utils library sources
set(AGPL_UTILS_SOURCES
    src/utils/allocator_helper.f90
    src/utils/console_colors.f90
    src/utils/logger.f90  # Fixed: Replaced '(XA)' with '(A)' format descriptor
    src/utils/operating_system.f90  # Fixed: Replaced SYSTEM calls with EXECUTE_COMMAND_LINE
    src/utils/random_numbers.f90
    src/utils/strings.f90
)

# Math library sources  
set(AGPL_MATH_SOURCES
    src/math/angles.f90
    src/math/coulomb_FG.f90
    src/math/double_factorial.f90
    src/math/fit_module.f90
    src/math/integration.f90
    src/math/Laguerre_polynomials.f90
    src/math/number_differences.f90
    src/math/spherical_Bessel_functions.f90
)

# Physics library sources
set(AGPL_PHYSICS_SOURCES
    src/physics/physical_constants.f90
    src/physics/quantum_numbers.f90
    src/physics/quantum_operators.f90
    src/physics/potential_partial_wave_caller.f90  # Re-enabled: EFT_pless dependency now available
    src/physics/scattering_NN.f90
    src/physics/potentials/av18.f90
    src/physics/potentials/av14.f90
    src/physics/potentials/EFT_pless.f90  # Re-enabled: operating_system dependency now fixed
)

# Create agpl_utils static library
add_library(agpl_utils STATIC ${AGPL_UTILS_SOURCES})
set_target_properties(agpl_utils PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    Fortran_MODULE_DIRECTORY "${FORTRAN_MODULE_DIRECTORY}"
)

# Create agpl_math static library
add_library(agpl_math STATIC ${AGPL_MATH_SOURCES})
target_link_libraries(agpl_math OpenMP::OpenMP_Fortran)
set_target_properties(agpl_math PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    Fortran_MODULE_DIRECTORY "${FORTRAN_MODULE_DIRECTORY}"
)

# Create agpl_physics static library
add_library(agpl_physics STATIC ${AGPL_PHYSICS_SOURCES})
# Physics depends on math (angles) and utils (allocator, quantum_numbers)
add_dependencies(agpl_physics agpl_math agpl_utils)
target_link_libraries(agpl_physics agpl_math agpl_utils)
set_target_properties(agpl_physics PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    Fortran_MODULE_DIRECTORY "${FORTRAN_MODULE_DIRECTORY}"
)

# Install targets
install(TARGETS agpl_utils agpl_math agpl_physics
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Install module files
install(DIRECTORY ${FORTRAN_MODULE_DIRECTORY}/
    DESTINATION include
    FILES_MATCHING PATTERN "*.mod"
)

# Optional: create a combined target for all libraries
add_custom_target(all_libraries 
    DEPENDS agpl_utils agpl_math agpl_physics
)

# Add examples subdirectory
add_subdirectory(examples)

# Print build information
message(STATUS "Building AGPL Fortran libraries:")
message(STATUS "  - agpl_utils: ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/")
message(STATUS "  - agpl_math: ${CMAKE_CURRENT_SOURCE_DIR}/src/math/")
message(STATUS "  - agpl_physics: ${CMAKE_CURRENT_SOURCE_DIR}/src/physics/")
message(STATUS "Output directory: ${CMAKE_BINARY_DIR}/lib/")
message(STATUS "Module directory: ${FORTRAN_MODULE_DIRECTORY}")