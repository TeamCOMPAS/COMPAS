
# use GNU C++ compiler by default
#
# can be overridden with CPP parameter
#
# e.g. make CPP=clang will use clang instead of g++
# (note uppercase 'CPP' and no whitespace around '=')

CPP := g++

# gsl directories
GSLINCDIR := /include
GSLLIBDIR := /lib

# boost directories
BOOSTINCDIR := /include
BOOSTLIBDIR := /lib

# hdf5 directories
HDF5INCDIR := /usr/include/hdf5/serial
HDF5LIBDIR := /usr/lib/x86_64-linux-gnu/hdf5/serial

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


CXXFLAGS := -std=c++11 -Wall -Woverloaded-virtual $(OPTFLAGS)
ICFLAGS := -I$(GSLINCDIR) -I$(BOOSTINCDIR) -I$(HDF5INCDIR) -I.

LIBS := -lm -lz -ldl -lpthread
GSLLIBS := -lgsl -lgslcblas
BOOSTLIBS := -lboost_filesystem -lboost_program_options -lboost_system
HDF5LIBS := -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
LFLAGS := -L$(GSLLIBDIR) -L$(BOOSTLIBDIR) -L$(HDF5LIBDIR) -Xlinker -rpath -Xlinker $(BOOSTLIBDIR) $(HDF5LIBS) $(LIBS) $(GSLLIBS) $(BOOSTLIBS)

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

%.o: %.cpp
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -c $?

.phony: clean static fast staticfast

fast: $(EXE)
staticfast:$(EXE)_STATIC

clean:
	rm -f $(OBJI) $(EXE) $(EXE)_STATIC
