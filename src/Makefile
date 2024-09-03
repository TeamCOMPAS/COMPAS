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
# Note: Adjust library paths using environment variables:
#   GSL_PREFIX, BOOST_PREFIX, HDF5_PREFIX
#   Example: make GSL_PREFIX=/custom/path
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

# Check if pkg-config is available and the libraries are registered
ifeq ($(shell which $(PKG_CONFIG) 2>/dev/null),)
  $(warning pkg-config not found, falling back to default paths)
  USE_PKG_CONFIG := 0
else
  ifeq ($(shell $(PKG_CONFIG) --exists gsl && echo yes),yes)
    USE_PKG_CONFIG := 1
  else
    $(warning gsl not found in pkg-config, falling back to default paths)
    USE_PKG_CONFIG := 0
  endif
endif

# Set include and lib directories based on pkg-config or prefix
ifeq ($(USE_PKG_CONFIG),1)
  GSL_CFLAGS := $(shell $(PKG_CONFIG) --cflags gsl)
  GSL_LIBS := $(shell $(PKG_CONFIG) --libs gsl)
else
  GSL_CFLAGS := -I$(GSL_PREFIX)/include
  GSL_LIBS := -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
endif

# Always use default paths for HDF5
HDF5_CFLAGS := -I$(HDF5_PREFIX)/include
HDF5_LIBS := -L$(HDF5_PREFIX)/lib -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

# Boost doesn't typically use pkg-config
BOOST_CFLAGS := -I$(BOOST_PREFIX)/include
BOOST_LIBS := -L$(BOOST_PREFIX)/lib -lboost_filesystem -lboost_program_options -lboost_system

EXE := COMPAS

# build COMPAS
ifeq ($(filter clean,$(MAKECMDGOALS)),)
  $(info Building $(EXE) with $(CPP))
endif

OPTFLAGS :=
ifneq ($(filter fast,$(MAKECMDGOALS)),)
  $(info Adding optimisation flags into the compilation - will take longer to build)
  OPTFLAGS += -march=native -O3
endif

ifneq ($(filter staticfast,$(MAKECMDGOALS)),)
  $(info Adding optimisation flags into the (static) compilation - will take longer to build)
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

.phony: clean static fast staticfast

fast: $(EXE)
staticfast: $(EXE)_STATIC

clean:
	rm -f $(OBJI) $(EXE) $(EXE)_STATIC
