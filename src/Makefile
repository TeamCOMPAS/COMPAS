
# use GNU C++ compiler by default
#
# can be overridden with CPP parameter
#
# e.g. make CPP=clang will use clang instaed of g++
# (note uppercase 'CPP' and no whitespace around '=')

CPP = g++
BOOST :=


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


CXXFLAGS := -g -std=c++11 -Wall $(OPTFLAGS)
ICFLAGS := -I$(BOOST)/include -I.
LFLAGS :=  -L$(BOOST)/lib -lgsl -lgslcblas -lstdc++ -lm -lz -ldl -lboost_filesystem -lboost_program_options -lboost_system -Xlinker -rpath -Xlinker $(BOOST)/lib
SOURCES :=						\
	profiling.cpp               \
	utils.cpp                   \
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

.cpp.o: $(SOURCES) $(INCL) Makefile
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -c $?

.phony: clean static fast staticfast

fast: $(EXE)
staticfast:$(EXE)_STATIC

clean:
	rm -f $(OBJI) $(EXE) $(EXE)_STATIC
