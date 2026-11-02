CPP = g++

UNAME_S := $(shell uname -s)
HOMEBREW_PREFIX := $(shell brew --prefix 2>/dev/null)

# Build configuration: DOCKER_BUILD=1 for Docker, otherwise local
ifeq ($(DOCKER_BUILD),1)
    # Docker-specific settings
    ODIR := ../obj
    BDIR := ../bin
    OPTFLAGS_SPECIFIC := -O3
    CXXFLAGS_SPECIFIC :=
    LFLAGS_SPECIFIC :=
else
    # Local-specific settings
    ODIR := obj
    BDIR := bin
    OPTFLAGS_SPECIFIC := -march=native -O3
    CXXFLAGS_SPECIFIC := -g -fnon-call-exceptions -Woverloaded-virtual
    LFLAGS_SPECIFIC := -rdynamic
endif

# Common variables for both builds
EXE := $(BDIR)/COMPAS

SOURCES :=						\
	profiling.cpp               \
	utils.cpp                   \
	yaml.cpp                   	\
	vector3d.cpp                \
								\
	Rand.cpp                    \
	Options.cpp                 \
	Log.cpp                     \
	Errors.cpp                  \
								\
	BaseStar.cpp                \
								\
	Star.cpp                    \
								\
	MainSequence.cpp            \
	MS_lte_07.cpp               \
	MS_gt_07.cpp                \
								\
	CH.cpp                      \
								\
	GiantBranch.cpp             \
	HG.cpp                      \
	FGB.cpp                     \
	CHeB.cpp                    \
	EAGB.cpp                    \
	TPAGB.cpp                   \
								\
	HeMS.cpp                    \
	HeHG.cpp                    \
	HeGB.cpp					\
								\
	Remnants.cpp				\
								\
	WhiteDwarfs.cpp				\
	HeWD.cpp					\
	COWD.cpp                    \
	ONeWD.cpp                   \
								\
	NS.cpp                      \
	BH.cpp                      \
	MR.cpp                      \
								\
	BinaryConstituentStar.cpp   \
	BaseBinaryStar.cpp          \
	CircumbinaryDiskFormation.cpp       \
	BinaryStar.cpp              \
								\
	main.cpp

# Optional spindler-c backend.  Default OFF keeps ordinary COMPAS builds
# independent of the external spindler-c and librinterpolate libraries.
# To enable the backend, either set this here or provide it on the make command
# line, e.g. make SPINDLER_C_ROOT=/path/to/spindler-c
SPINDLER_C_ROOT ?=

ifneq ($(strip $(SPINDLER_C_ROOT)),)
SOURCES += CircumbinaryDiskEvolutionSpindler.cpp
else
SOURCES += CircumbinaryDiskEvolutionUnavailable.cpp
endif

OBJS := $(patsubst %.cpp,$(ODIR)/%.o,$(SOURCES))

# GSL, Boost, and HDF5 directories
GSLINCDIR ?= /include
GSLLIBDIR ?= /lib

# boost directories
BOOSTINCDIR ?= /include
BOOSTLIBDIR ?= /lib

# hdf5 directories
HDF5INCDIR ?= /usr/include/hdf5/serial
HDF5LIBDIR ?= /usr/lib/$(shell $(CPP) -print-multiarch)/hdf5/serial

# spindler-c include/link flags are populated only when SPINDLER_C_ROOT is
# nonempty.
SPINDLER_C_CPPFLAGS :=
SPINDLER_C_LDFLAGS :=
SPINDLER_C_LIBS :=

ifneq ($(strip $(SPINDLER_C_ROOT)),)

SPINDLER_C_INCDIR ?= $(SPINDLER_C_ROOT)/include
SPINDLER_C_LIBDIR ?= $(SPINDLER_C_ROOT)/lib

# librinterpolate is a dependency of spindler-c, linked as -lrinterpolate.
# Leave RINTERPOLATE_ROOT empty if librinterpolate is installed in a standard
# compiler/linker search path. Otherwise, set it here or on the make command
# line, e.g. make SPINDLER_C_ROOT=/path/to/spindler-c RINTERPOLATE_ROOT=/path/to/rinterpolate
RINTERPOLATE_ROOT ?=
RINTERPOLATE_INCDIR ?=
RINTERPOLATE_LIBDIR ?=

ifneq ($(strip $(RINTERPOLATE_ROOT)),)
ifeq ($(strip $(RINTERPOLATE_INCDIR)),)
RINTERPOLATE_INCDIR := $(RINTERPOLATE_ROOT)/include
endif
ifeq ($(strip $(RINTERPOLATE_LIBDIR)),)
RINTERPOLATE_LIBDIR := $(RINTERPOLATE_ROOT)/lib
endif
endif

SPINDLER_C_CPPFLAGS += -DCOMPAS_ENABLE_SPINDLER_C
SPINDLER_C_CPPFLAGS += -I$(SPINDLER_C_INCDIR)
ifneq ($(strip $(RINTERPOLATE_INCDIR)),)
SPINDLER_C_CPPFLAGS += -I$(RINTERPOLATE_INCDIR)
endif
SPINDLER_C_CPPFLAGS += -DSPINDLER_C_ROOT=\"$(SPINDLER_C_ROOT)\"
SPINDLER_C_LDFLAGS += -L$(SPINDLER_C_LIBDIR) \
                      -Wl,-rpath,$(SPINDLER_C_LIBDIR)
ifneq ($(strip $(RINTERPOLATE_LIBDIR)),)
SPINDLER_C_LDFLAGS += -L$(RINTERPOLATE_LIBDIR) \
                      -Wl,-rpath,$(RINTERPOLATE_LIBDIR)
endif
SPINDLER_C_LIBS += -lspindler -lrinterpolate

endif

# Homebrew defaults for macOS local builds
ifeq ($(UNAME_S),Darwin)
ifneq ($(wildcard $(HOMEBREW_PREFIX)/opt/boost/include/boost/version.hpp),)
    GSLINCDIR := $(HOMEBREW_PREFIX)/opt/gsl/include
    GSLLIBDIR := $(HOMEBREW_PREFIX)/opt/gsl/lib
    BOOSTINCDIR := $(HOMEBREW_PREFIX)/opt/boost/include
    BOOSTLIBDIR := $(HOMEBREW_PREFIX)/opt/boost/lib
    HDF5INCDIR := $(HOMEBREW_PREFIX)/opt/hdf5/include
    HDF5LIBDIR := $(HOMEBREW_PREFIX)/opt/hdf5/lib
endif
endif

# Define flags
OPTFLAGS :=
ifneq ($(filter fast staticfast,$(MAKECMDGOALS)),)
    $(info Adding optimisation flags into the compilation - will take longer to build)
    OPTFLAGS += $(OPTFLAGS_SPECIFIC)
endif

# might need -Wno-vla-cxx-extension for clang>18?
EXTRA_WARN_SUPPRESS := -Wno-vla-extension 

CXXFLAGS := -std=c++17 -Wall $(OPTFLAGS) $(CXXFLAGS_SPECIFIC) $(EXTRA_WARN_SUPPRESS) $(SPINDLER_C_CPPFLAGS)
ICFLAGS := -I$(GSLINCDIR) -I$(BOOSTINCDIR) -I$(HDF5INCDIR) -I.

LIBS := -lm -lz -ldl -lpthread
GSLLIBS := -lgsl -lgslcblas
BOOSTLIBS := -lboost_filesystem -lboost_program_options -lboost_system
ifeq ($(UNAME_S),Darwin)
ifeq ($(wildcard $(BOOSTLIBDIR)/libboost_system*),)
    BOOSTLIBS := -lboost_filesystem -lboost_program_options
endif
endif
HDF5LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
LFLAGS := -L$(GSLLIBDIR) -L$(BOOSTLIBDIR) -L$(HDF5LIBDIR) $(SPINDLER_C_LDFLAGS) -Xlinker -rpath -Xlinker $(BOOSTLIBDIR) $(HDF5LIBS) $(LIBS) $(GSLLIBS) $(BOOSTLIBS) $(SPINDLER_C_LIBS) $(LFLAGS_SPECIFIC)

.PHONY: all static fast staticfast clean link

# Main targets
all: $(EXE) link

static: $(EXE)_STATIC

fast: all

staticfast: static

$(EXE): $(OBJS)
	@mkdir -p $(BDIR)
	$(CPP) $(OBJS) $(LFLAGS) -o $@

$(EXE)_STATIC: $(OBJS)
	@mkdir -p $(BDIR)
	$(CPP) $(OBJS) $(LFLAGS) -static -o $@

# Create a symbolic link for backward compatibility
link: $(EXE)
	@ln -sf $(EXE) COMPAS

# Pattern rule to compile source files into object files
$(ODIR)/%.o: %.cpp
	@mkdir -p $(ODIR)
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -o $@ -c $<

# Clean-up rule
clean:
	@echo "Removing generated files..."
	@rm -rf $(ODIR) $(BDIR) COMPAS
