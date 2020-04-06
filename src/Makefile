
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

CXXFLAGS := -g -std=c++11 -Wall -march=native -O3
ICFLAGS := -I$(BOOST)/include -I.
LFLAGS :=  -L$(BOOST)/lib -lgsl -lgslcblas -lstdc++ -lm -lz -ldl -lboost_filesystem -lboost_program_options -lboost_system -Xlinker -rpath -Xlinker $(BOOST)/lib
SOURCES :=                                                               \
			utils.cpp					\
									\
			Options.cpp					\
			Log.cpp						\
			Rand.cpp					\
			Errors.cpp					\
									\
			BaseStar.cpp					\
									\
			Star.cpp					\
									\
			MainSequence.cpp				\
			MS_lte_07.cpp					\
			MS_gt_07.cpp					\
									\
			CH.cpp						\
									\
			GiantBranch.cpp					\
			HG.cpp						\
			FGB.cpp						\
			CHeB.cpp					\
			EAGB.cpp					\
			TPAGB.cpp					\
									\
			HeMS.cpp					\
			HeHG.cpp					\
			HeGB.cpp					\
									\
			HeWD.cpp					\
			COWD.cpp					\
			ONeWD.cpp					\
									\
			NS.cpp						\
			BH.cpp						\
			MR.cpp						\
									\
			AIS.cpp						\
									\
			BinaryConstituentStar.cpp			\
			BaseBinaryStar.cpp				\
			BinaryStar.cpp					\
									\
			main.cpp

OBJI := $(SOURCES:.cpp=.o)


all: $(EXE)
	@echo $(OBJI)

.cpp.o: $(SOURCES)
	$(CPP) $(CXXFLAGS) $(ICFLAGS) -c $?


$(EXE): $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) $(OBJI) $(LFLAGS) -o $@

.phony: clean static
clean:
	rm -f $(OBJI)

static: $(EXE)_STATIC
	@echo $(OBJI)


$(EXE)_STATIC: $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) $(OBJI) $(LFLAGS) -static -o $@
