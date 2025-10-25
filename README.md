# AGPL Fortran Libraries Build System

## Overview

This CMakeLists.txt file creates three static libraries from your Fortran source code:

- **libagpl_utils.a** - Utility functions (strings, colors, allocators, random numbers)
- **libagpl_math.a** - Mathematical functions (angles, polynomials, special functions)
- **libagpl_physics.a** - Physics calculations (constants, quantum numbers, potentials)

## Successfully Built Modules

### agpl_utils library:
- allocator_helper.f90 → reallocate_utils module
- console_colors.f90 → console_colors module  
- logger.f90 → log module ✅ **Fixed**
- operating_system.f90 → operating_system_linux module ✅ **Fixed**
- random_numbers.f90 → random_numbers module
- strings.f90 → strings_utils module

### agpl_math library:
- angles.f90 → angles module
- coulomb_FG.f90 → gsl_coulomb module
- double_factorial.f90 (provides double factorial functions)
- fit_module.f90 → fit_module module
- Laguerre_polynomials.f90 → laguerre_polynomial_mod module
- number_differences.f90 → number_differences module
- spherical_Bessel_functions.f90 → gsl_bessel module

### agpl_physics library:
- physical_constants.f90 → physical_constants module
- quantum_numbers.f90 → quantum_numbers module
- quantum_operators.f90 → quantum_operators module
- potential_partial_wave_caller.f90 → potentials module ✅ **Fixed**
- scattering_NN.f90 → scattering module ✅ **Fixed**
- potentials/av18.f90 → av18 module
- potentials/av14.f90 → av14 module
- potentials/EFT_pless.f90 → eft_pless module ✅ **Fixed**

## Build Instructions

```bash
# Create build directory
mkdir -p build && cd build

# Configure with gfortran (automatically detected)
cmake ..

# Or using the Makefile wrapper (recommended)
make configure

# Build all libraries
make -j$(nproc)

# Or using the wrapper
make compile
```

## Makefile Wrapper

The project includes a convenient Makefile wrapper that provides easy-to-use commands:

```bash
make help              # Show all available commands
make all               # Configure and build everything
make configure         # Run CMake configuration
make compile           # Build all libraries
make utils/math/physics # Build individual libraries
make clean             # Clean build artifacts
make info              # Show build status
make compile-commands  # Show compilation commands database
```

## Compile Commands Database

The build system automatically generates a `compile_commands.json` file containing all compilation commands. This is useful for:

- **Language servers**: For IDEs and editors like VS Code, Vim, Emacs
- **Static analysis tools**: For tools like clang-tidy, cppcheck
- **Build debugging**: To see exact compiler flags and commands used

Access the compile commands with:
```bash
make compile-commands  # Pretty-printed view
cat build/compile_commands.json  # Raw JSON
```

## Library Files Location

- Static libraries: `build/lib/libagpl_*.a`
- Module files: `build/modules/*.mod`

## Usage Example

To use these libraries in your projects:

```bash
# Compile your program linking against the libraries
gfortran your_program.f90 -I/path/to/modules -L/path/to/lib -lagpl_physics -lagpl_math -lagpl_utils -fopenmp
```

## Files Currently Excluded (Need Fixes)

### ~~src/utils/logger.f90~~ ✅ **FIXED**
**Issue**: X descriptor format needs leading space count  
**Fix**: ✅ **COMPLETED** - Changed `'(XA)'` to `'(A)'` on lines 340 and 342

### ~~src/utils/operating_system.f90~~ ✅ **FIXED**
**Issue**: SYSTEM intrinsic call syntax  
**Fix**: ✅ **COMPLETED** - Replaced all `CALL SYSTEM(command)` with `CALL EXECUTE_COMMAND_LINE(command)` and `CALL SYSTEM(command, STATUS=ios)` with `CALL EXECUTE_COMMAND_LINE(command, EXITSTAT=ios)`

### src/math/integration.f90
**Issue**: DFLOAT intrinsic not recognized
**Fix**: Replace `DFLOAT(N)` with `REAL(N, KIND=REAL64)` on line 248

### ~~src/physics/scattering_NN.f90~~ ✅ **FIXED**
**Issue**: DOUBLE COMPLEX is non-standard  
**Fix**: ✅ **COMPLETED** - Replaced `DOUBLE COMPLEX` with `COMPLEX(KIND=REAL64)`, `CDEXP` with `EXP`, `CDSQRT` with `SQRT`, `DREAL` with `REAL`, `DIMAG` with `AIMAG`

### ~~src/physics/potentials/EFT_pless.f90~~ ✅ **FIXED**
**Issue**: Depends on operating_system module which had compilation issues  
**Fix**: ✅ **COMPLETED** - operating_system.f90 is now fixed, EFT_pless.f90 successfully re-enabled

### ~~src/physics/potential_partial_wave_caller.f90~~ ✅ **FIXED**
**Issue**: Depends on EFT_pless module  
**Fix**: ✅ **COMPLETED** - EFT_pless.f90 is now fixed, potential_partial_wave_caller.f90 successfully re-enabled

## Dependencies

- **OpenMP**: Required for math library (parallel integration)
- **gfortran**: GNU Fortran compiler (tested with version 14.2.0)
- **CMake**: Version 3.12 or higher

## Compiler Flags

- **Standard**: Fortran 2008 (`-std=f2008`)
- **Position Independent Code**: `-fPIC`
- **Release Mode**: `-O3 -march=native`
- **Debug Mode**: `-g -O0 -fcheck=all -fbacktrace -Wall`

## Installation

The CMake configuration supports installation:

```bash
make install
```

This will install:
- Libraries to `CMAKE_INSTALL_PREFIX/lib/`
- Module files to `CMAKE_INSTALL_PREFIX/include/`

## Module Dependencies

The libraries have been organized with proper dependency handling:
- `agpl_physics` depends on `agpl_math` (for angles module)
- `agpl_physics` depends on `agpl_utils` (for allocator module)
- `agpl_math` has OpenMP dependency for parallel processing

## Notes

- All successfully compiled modules are available in `build/modules/`
- The build system handles module dependencies automatically
- Libraries are built as static archives (`.a` files)
- Position-independent code is enabled for potential shared library use

To include all files (after fixing the issues), uncomment the relevant lines in CMakeLists.txt and rebuild.