cmake_minimum_required(VERSION 3.12)

# Set the example output directory
set(EXAMPLE_OUTPUT_DIR ${CMAKE_BINARY_DIR}/examples)

# Create the examples output directory
file(MAKE_DIRECTORY ${EXAMPLE_OUTPUT_DIR})

# Include the modules directory for Fortran modules
include_directories(${FORTRAN_MODULE_DIRECTORY})

# Example: Operating System Module Demo
add_executable(example_operating_system
    example_operating_system.f90
)

# Link against the agpl_utils library (which contains operating_system module)
target_link_libraries(example_operating_system agpl_utils)

# Set output directory for examples
set_target_properties(example_operating_system PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_DIR}
    Fortran_MODULE_DIRECTORY ${FORTRAN_MODULE_DIRECTORY}
)

# Add dependency to ensure libraries are built first
add_dependencies(example_operating_system agpl_utils)

# Print information about examples
message(STATUS "Examples will be built in: ${EXAMPLE_OUTPUT_DIR}")

# Optional: Create a target to build all examples
add_custom_target(examples
    DEPENDS example_operating_system
    COMMENT "Building all examples"
)

# Install examples (optional)
install(TARGETS example_operating_system
    RUNTIME DESTINATION examples
)