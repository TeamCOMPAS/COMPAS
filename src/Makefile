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
	BinaryStar.cpp              \
								\
	main.cpp

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

# might need -Wno-vla-cxx-extension for clang>18?
EXTRA_WARN_SUPPRESS := -Wno-vla-extension 

CXXFLAGS := -std=c++17 -Wall $(OPTFLAGS) $(CXXFLAGS_SPECIFIC) $(EXTRA_WARN_SUPPRESS)
ICFLAGS := -I$(GSLINCDIR) -I$(BOOSTINCDIR) -I$(HDF5INCDIR) -I.

LIBS := -lm -lz -ldl -lpthread
GSLLIBS := -lgsl -lgslcblas
BOOSTLIBS := -lboost_filesystem -lboost_program_options -lboost_system
HDF5LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
LFLAGS := -L$(GSLLIBDIR) -L$(BOOSTLIBDIR) -L$(HDF5LIBDIR) -Xlinker -rpath -Xlinker $(BOOSTLIBDIR) $(HDF5LIBS) $(LIBS) $(GSLLIBS) $(BOOSTLIBS) $(LFLAGS_SPECIFIC)

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
