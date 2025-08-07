
# use GNU C++ compiler by default
#
# can be overridden with CPP parameter
#
# e.g. make CPP=clang will use clang instead of g++
# (note uppercase 'CPP' and no whitespace around '=')

CPP := g++

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
SOURCES := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,$(ODIR)/%.o,$(SOURCES))

# GSL, Boost, and HDF5 directories
GSLINCDIR := /include
GSLLIBDIR := /lib

# boost directories
BOOSTINCDIR := /include
BOOSTLIBDIR := /lib

# hdf5 directories
HDF5INCDIR := /usr/include/hdf5/serial
HDF5LIBDIR := /usr/lib/x86_64-linux-gnu/hdf5/serial

# Define flags
OPTFLAGS :=
ifneq ($(filter fast staticfast,$(MAKECMDGOALS)),)
    $(info Adding optimisation flags into the compilation - will take longer to build)
    OPTFLAGS += $(OPTFLAGS_SPECIFIC)
endif

CXXFLAGS := -std=c++17 -Wall $(OPTFLAGS) $(CXXFLAGS_SPECIFIC)
ICFLAGS := -I$(GSLINCDIR) -I$(BOOSTINCDIR) -I$(HDF5INCDIR) -I.

LIBS := -lm -lz -ldl -lpthread
GSLLIBS := -lgsl -lgslcblas
BOOSTLIBS := -lboost_filesystem -lboost_program_options -lboost_system
HDF5LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
LFLAGS := -L$(GSLLIBDIR) -L$(BOOSTLIBDIR) -L$(HDF5LIBDIR) -Xlinker -rpath -Xlinker $(BOOSTLIBDIR) $(HDF5LIBS) $(LIBS) $(GSLLIBS) $(BOOSTLIBS) $(LFLAGS_SPECIFIC)

.PHONY: all static fast staticfast clean

# Main targets
all: $(EXE)

static: $(EXE)_STATIC

fast: all

staticfast: static

$(EXE): $(OBJS)
	@mkdir -p $(BDIR)
	$(CPP) $(OBJS) $(LFLAGS) -o $@

$(EXE)_STATIC: $(OBJS)
	@mkdir -p $(BDIR)
	$(CPP) $(OBJS) $(LFLAGS) -static -o $@

# Pattern rule to compile source files into object files
$(ODIR)/%.o: %.cpp
	@mkdir -p $(ODIR)
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -o $@ -c $<

# Clean-up rule
clean:
	@echo "Removing generated files..."
	@rm -rf $(ODIR) $(BDIR)