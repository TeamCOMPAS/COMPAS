#Use GNU C++ compiler
CPP = g++

#Check compiler version - need GCC/g++ 7.2.0 or greater
CPPVERSION_1_EQ_7 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \== 7)
CPPVERSION_2_GTEQ_2 := $(shell expr `gcc -dumpversion | cut -f2 -d.` \>= 2)
CPPVERSION_1_GTEQ_7 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 7)

OK = 0
ifeq "$(CPPVERSION_1_EQ_7)" "1"
	ifeq "$(CPPVERSION_2_GTEQ_2)" "1"
		OK = 1
	endif
else
	ifeq "$(CPPVERSION_1_GTEQ_7)" "1"
		OK = 1
	endif
endif

#if compiler version is good, compile
ifeq "$(OK)" "1"


CFLAGS = -g


BOOST =
CFLAGS =  -I.
FCFLAGS =
ICFLAGS = -I$(BOOST)/include -std=gnu++11
LFLAGS =  -L$(BOOST)/lib -lgsl -lgslcblas -lm -lz -ldl -march=k8 -lboost_filesystem -lboost_program_options -lboost_system -std=gnu++11 -O3
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

else

#bad compiler version - announce error
all:
	@echo "Need GCC/g++ version 7.2.0 or greater to compile COMPAS"

endif

