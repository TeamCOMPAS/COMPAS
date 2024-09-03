# Makefile Usage Guide
#
# Build the project:
#   make            - Compiles the project using the default compiler (g++).
#
# Use a different compiler:
#   make CPP=clang  - Compiles using clang instead of g++.
#
# Enable optimization:
#   make fast       - Compiles with optimization flags for faster execution.
#
# Build a static executable:
#   make static     - Builds a static version of the executable.
#
# Clean the build:
#   make clean      - Removes all compiled objects and executables.
#
# Parallel build:
#   make -j $(nproc) -f Makefile
#                  - Compiles the project using all available CPU cores for faster build times.
#                    '-j $(nproc)' specifies the number of jobs to run simultaneously.
#
# Library Configuration:
#   Adjust library paths using environment variables:
#   GSL_PREFIX, BOOST_PREFIX, HDF5_PREFIX
#   Example: make GSL_PREFIX=/custom/path
#
# HDF5 Configuration:
#   Ensure HDF5_PREFIX is set correctly:
#   make HDF5_PREFIX=/path/to/hdf5
#   or set HDF5_DIR environment variable before running make.
#
# Use pkg-config for dependencies if available; otherwise, default paths are used.

# use GNU C++ compiler by default
CPP := g++

# Use environment variables for library paths, with defaults
GSL_PREFIX ?= /usr
BOOST_PREFIX ?= /usr
HDF5_PREFIX ?= /usr

# Use pkg-config if available
PKG_CONFIG ?= pkg-config

# Check if pkg-config is available
ifeq ($(shell which $(PKG_CONFIG) 2>/dev/null),)
  $(warning pkg-config not found, falling back to default paths)
  USE_PKG_CONFIG := 0
else
  USE_PKG_CONFIG := 1
endif

# Function to set flags using pkg-config or defaults
define set_pkg_config_flags
  ifeq ($(USE_PKG_CONFIG),1)
    ifeq ($(shell $(PKG_CONFIG) --exists $(1) && echo yes),yes)
      $(2)_CFLAGS := $(shell $(PKG_CONFIG) --cflags $(1))
      $(2)_LIBS := $(shell $(PKG_CONFIG) --libs $(1))
    else
      $(warning $(1) not found in pkg-config, falling back to default paths)
      $(2)_CFLAGS := -I$($(2)_PREFIX)/include
      $(2)_LIBS := -L$($(2)_PREFIX)/lib $(3)
    endif
  else
    $(2)_CFLAGS := -I$($(2)_PREFIX)/include
    $(2)_LIBS := -L$($(2)_PREFIX)/lib $(3)
  endif
endef

# Use pkg-config for GSL, Boost, and HDF5
$(eval $(call set_pkg_config_flags,gsl,GSL,-lgsl -lgslcblas))
$(eval $(call set_pkg_config_flags,boost,BOOST,-lboost_filesystem -lboost_program_options -lboost_system))
$(eval $(call set_pkg_config_flags,hdf5,HDF5,-lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5))

EXE := COMPAS

# Build COMPAS
ifeq ($(filter clean,$(MAKECMDGOALS)),)
  $(info Building $(EXE) with $(CPP))
endif

OPTFLAGS :=
ifneq ($(filter fast,$(MAKECMDGOALS)),)
  $(info Adding optimization flags into the compilation - will take longer to build)
  OPTFLAGS += -march=native -O3
endif

ifneq ($(filter staticfast,$(MAKECMDGOALS)),)
  $(info Adding optimization flags into the (static) compilation - will take longer to build)
  OPTFLAGS += -march=native -O3
endif

CXXFLAGS := -std=c++11 -Wall $(OPTFLAGS)
ICFLAGS := $(GSL_CFLAGS) $(BOOST_CFLAGS) $(HDF5_CFLAGS) -I.

LIBS := -lm -lz -ldl -lpthread
LFLAGS := $(GSL_LIBS) $(BOOST_LIBS) $(HDF5_LIBS) $(LIBS)

SOURCES := profiling.cpp utils.cpp yaml.cpp vector3d.cpp \
           Rand.cpp Options.cpp Log.cpp Errors.cpp \
           BaseStar.cpp Star.cpp \
           MainSequence.cpp MS_lte_07.cpp MS_gt_07.cpp \
           CH.cpp \
           GiantBranch.cpp HG.cpp FGB.cpp CHeB.cpp EAGB.cpp TPAGB.cpp \
           HeMS.cpp HeHG.cpp HeGB.cpp \
           Remnants.cpp \
           WhiteDwarfs.cpp HeWD.cpp COWD.cpp ONeWD.cpp \
           NS.cpp BH.cpp MR.cpp \
           BinaryConstituentStar.cpp BaseBinaryStar.cpp BinaryStar.cpp \
           main.cpp

OBJI := $(SOURCES:.cpp=.o)

# Create the list of header files, and remove
# main.h from this auto-generated list
INCL := $(SOURCES:.cpp=.h)
INCL := $(filter-out main.h,$(INCL))

all: $(EXE)
	@echo $(OBJI)

$(EXE): $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) $(OBJI) $(LFLAGS) -o $@

static: $(EXE)_STATIC
	@echo $(OBJI)

$(EXE)_STATIC: $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) $(OBJI) $(LFLAGS) -static -o $@

.cpp.o: $(SOURCES) $(INCL) Makefile
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -c $?

.PHONY: clean static fast staticfast

fast: $(EXE)
staticfast: $(EXE)_STATIC

clean:
	rm -f $(OBJI) $(EXE) $(EXE)_STATIC

# Debug target to print variables
print-%:
	@echo $* = $($*)