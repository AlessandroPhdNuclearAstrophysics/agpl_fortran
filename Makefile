# Makefile wrapper for CMake-based AGPL Fortran Libraries
# This provides convenient make targets that use CMake underneath

# Default build type
BUILD_TYPE ?= Release

# Build directory
BUILD_DIR = build

# Number of parallel jobs
JOBS ?= $(shell nproc)

# Compiler selection
FC := $(shell which gfortran 2>/dev/null || which ifort 2>/dev/null || which ifx 2>/dev/null || echo gfortran)

# Default target
.PHONY: all
all: configure compile

# Help target
.PHONY: help
help:
	@echo "AGPL Fortran Libraries - Makefile Wrapper"
	@echo ""
	@echo "Available targets:"
	@echo "  all          - Configure and build all libraries (default)"
	@echo "  configure    - Run CMake configuration"
	@echo "  compile      - Build all libraries"
	@echo "  clean        - Remove build artifacts"
	@echo "  distclean    - Remove entire build directory"
	@echo "  install      - Install libraries and modules"
	@echo "  test         - Run tests (if available)"
	@echo "  libs         - Build only the libraries"
	@echo "  utils        - Build only agpl_utils library"
	@echo "  math         - Build only agpl_math library"
	@echo "  physics      - Build only agpl_physics library"
	@echo "  examples     - Build example programs"
	@echo "  debug        - Build with debug flags"
	@echo "  release      - Build with release optimization"
	@echo "  info         - Show build information"
	@echo "  list-modules - List generated module files"
	@echo "  list-libs    - List generated library files"
	@echo "  compile-commands - Show compilation commands database"
	@echo ""
	@echo "Variables:"
	@echo "  BUILD_TYPE   - Release (default) or Debug"
	@echo "  FC           - Fortran compiler (default: gfortran)"
	@echo "  JOBS         - Number of parallel jobs (default: $(JOBS))"
	@echo ""
	@echo "Examples:"
	@echo "  make                    # Configure and build with default settings"
	@echo "  make debug              # Build in debug mode"
	@echo "  make FC=ifort           # Use Intel Fortran compiler"
	@echo "  make BUILD_TYPE=Debug   # Explicit debug build"
	@echo "  make clean compile      # Clean and rebuild"

# Create build directory
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Configure with CMake
.PHONY: configure
configure: $(BUILD_DIR)
	@echo "Configuring with CMake..."
	@echo "  Fortran Compiler: $(FC)"
	@echo "  Build Type: $(BUILD_TYPE)"
	CMAKE_Fortran_COMPILER=$(FC) FC=$(FC) cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_Fortran_COMPILER=$(FC)

# Build all targets
.PHONY: compile
compile:
	@if [ ! -f $(BUILD_DIR)/Makefile ]; then \
		echo "Build not configured. Running configure first..."; \
		$(MAKE) configure; \
	fi
	@echo "Building libraries with $(JOBS) parallel jobs..."
	$(MAKE) -C $(BUILD_DIR) -j$(JOBS)

# Individual library targets
.PHONY: libs utils math physics
libs: compile

utils:
	@if [ ! -f $(BUILD_DIR)/Makefile ]; then \
		echo "Build not configured. Running configure first..."; \
		$(MAKE) configure; \
	fi
	@echo "Building agpl_utils library..."
	$(MAKE) -C $(BUILD_DIR) agpl_utils

math:
	@if [ ! -f $(BUILD_DIR)/Makefile ]; then \
		echo "Build not configured. Running configure first..."; \
		$(MAKE) configure; \
	fi
	@echo "Building agpl_math library..."
	$(MAKE) -C $(BUILD_DIR) agpl_math

physics:
	@if [ ! -f $(BUILD_DIR)/Makefile ]; then \
		echo "Build not configured. Running configure first..."; \
		$(MAKE) configure; \
	fi
	@echo "Building agpl_physics library..."
	$(MAKE) -C $(BUILD_DIR) agpl_physics

examples:
	@if [ ! -f $(BUILD_DIR)/Makefile ]; then \
		echo "Build not configured. Running configure first..."; \
		$(MAKE) configure; \
	fi
	@echo "Building example programs..."
	$(MAKE) -C $(BUILD_DIR) example_operating_system

# Build type shortcuts
.PHONY: debug release
debug:
	@$(MAKE) BUILD_TYPE=Debug all

release:
	@$(MAKE) BUILD_TYPE=Release all

# Clean targets
.PHONY: clean distclean
clean:
	@if [ -d $(BUILD_DIR) ]; then \
		echo "Cleaning build artifacts..."; \
		$(MAKE) -C $(BUILD_DIR) clean; \
	else \
		echo "No build directory found."; \
	fi

distclean:
	@echo "Removing entire build directory..."
	rm -rf $(BUILD_DIR)

# Install target
.PHONY: install
install: compile
	@echo "Installing libraries and modules..."
	$(MAKE) -C $(BUILD_DIR) install

# Test target (placeholder)
.PHONY: test
test: compile
	@echo "Running tests..."
	@if [ -f $(BUILD_DIR)/CTestTestfile.cmake ]; then \
		ctest --test-dir $(BUILD_DIR); \
	else \
		echo "No tests configured."; \
	fi

# Information targets
.PHONY: info list-modules list-libs
info:
	@echo "AGPL Fortran Libraries Build Information"
	@echo "========================================"
	@echo "Fortran Compiler: $(shell which gfortran 2>/dev/null || echo $(FC))"
	@echo "Build Type: $(BUILD_TYPE)"
	@echo "Build Directory: $(BUILD_DIR)"
	@echo "Parallel Jobs: $(JOBS)"
	@echo ""
	@if [ -d $(BUILD_DIR) ]; then \
		echo "Build Status: Configured"; \
		if [ -f $(BUILD_DIR)/lib/libagpl_utils.a ]; then echo "  ✓ agpl_utils library built"; fi; \
		if [ -f $(BUILD_DIR)/lib/libagpl_math.a ]; then echo "  ✓ agpl_math library built"; fi; \
		if [ -f $(BUILD_DIR)/lib/libagpl_physics.a ]; then echo "  ✓ agpl_physics library built"; fi; \
	else \
		echo "Build Status: Not configured"; \
	fi

list-modules:
	@echo "Generated Fortran Module Files:"
	@echo "==============================="
	@if [ -d $(BUILD_DIR)/modules ]; then \
		ls -la $(BUILD_DIR)/modules/*.mod 2>/dev/null || echo "No module files found."; \
	else \
		echo "No modules directory found. Build first."; \
	fi

list-libs:
	@echo "Generated Library Files:"
	@echo "======================="
	@if [ -d $(BUILD_DIR)/lib ]; then \
		ls -la $(BUILD_DIR)/lib/*.a 2>/dev/null || echo "No library files found."; \
	else \
		echo "No lib directory found. Build first."; \
	fi

# Rebuild target
.PHONY: rebuild
rebuild: clean compile

# Quick development targets
.PHONY: quick dev
quick: configure
	$(MAKE) -C $(BUILD_DIR) -j$(JOBS)

dev: debug

# Documentation target
.PHONY: docs
docs:
	@echo "Documentation files:"
	@echo "  README_BUILD.md - Build instructions and troubleshooting"
	@echo "  CMakeLists.txt  - CMake configuration"
	@echo "  Makefile        - This makefile wrapper"
	@if [ -f $(BUILD_DIR)/compile_commands.json ]; then \
		echo "  $(BUILD_DIR)/compile_commands.json - Compilation commands database"; \
	fi

# Show compile commands
.PHONY: compile-commands show-commands
compile-commands show-commands:
	@if [ -f $(BUILD_DIR)/compile_commands.json ]; then \
		echo "Compile commands database:"; \
		echo "========================"; \
		cat $(BUILD_DIR)/compile_commands.json | jq -r '.[] | "File: " + .file + "\nCommand: " + .command + "\n"' 2>/dev/null || \
		python3 -m json.tool $(BUILD_DIR)/compile_commands.json 2>/dev/null || \
		cat $(BUILD_DIR)/compile_commands.json; \
	else \
		echo "No compile commands database found."; \
		echo "Run 'make configure' or 'make all' first."; \
	fi

# Force reconfigure
.PHONY: reconfigure
reconfigure:
	@echo "Forcing reconfiguration..."
	rm -f $(BUILD_DIR)/CMakeCache.txt
	$(MAKE) configure

# Show compiler version
.PHONY: compiler-info
compiler-info:
	@echo "Fortran Compiler Information:"
	@$(FC) --version 2>/dev/null || echo "Compiler $(FC) not found"
	@echo ""
	@echo "Available Fortran compilers:"
	@which gfortran 2>/dev/null && echo "  ✓ gfortran" || echo "  ✗ gfortran"
	@which ifort 2>/dev/null && echo "  ✓ ifort" || echo "  ✗ ifort"
	@which ifx 2>/dev/null && echo "  ✓ ifx" || echo "  ✗ ifx"