
# use GNU C++ compiler by default
#
# can be overridden with CPP parameter
#
# e.g. make CPP=clang will use clang instaed of g++
# (note uppercase 'CPP' and no whitespace around '=')

CPP = g++
cpp = $(CPP)


# build COMPAS

$(info Building COMPAS with $(CPP))

CFLAGS = -g


BOOST =
CFLAGS =  -I.
FCFLAGS =
ICFLAGS = -I$(BOOST)/include -Wall -std=c++11
LFLAGS =  -L$(BOOST)/lib -lgsl -lgslcblas -lstdc++ -lm -lz -ldl -march=k8 -lboost_filesystem -lboost_program_options -lboost_system -std=c++11 -O3
SOURCES = 								\
			utils.cpp					\
										\
			Options.cpp					\
			Log.cpp						\
			Rand.cpp					\
			Errors.cpp					\
										\
			BaseStar.cpp				\
										\
			Star.cpp					\
										\
			MainSequence.cpp			\
			MS_lte_07.cpp				\
			MS_gt_07.cpp				\
										\
			CH.cpp						\
										\
			GiantBranch.cpp				\
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
			BinaryConstituentStar.cpp	\
			BaseBinaryStar.cpp			\
			BinaryStar.cpp				\
										\
			main.cpp

OBJI = $(SOURCES:.cpp=.o)


all: COMPAS
	@echo $(OBJI)

.cpp.o:
	$(CPP) -c $(CFLAGS) $(FCFLAGS) $(ICFLAGS) $?


COMPAS: $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) -o COMPAS $(OBJI) $(ICFLAGS) $(LFLAGS)

clean:
	rm -f *.o

static: COMPAS_STATIC
	@echo $OBJI

COMPAS_STATIC: $(OBJI)
	@echo $(SOURCES)
	@echo $(OBJI)
	$(CPP) -o COMPAS $(OBJI) $(ICFLAGS) $(LFLAGS) -static


