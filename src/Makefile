
# use GNU C++ compiler by default
#
# can be overridden with CPP parameter
#
# e.g. make CPP=clang will use clang instaed of g++
# (note uppercase 'CPP' and no whitespace around '=')

CPP = g++
cpp = $(CPP)


# we check what compiler is being used, and the version
# if the compiler and/or version is not what we expect,
# COMPAS is not built
#
# default is to check - can be overridden with CHECK parameter
#
# e.g. make CHECK=off (or make CHECK=OFF) to disable the check
# (note uppercase 'CHECK' and no whitespace around '=')

CHECK = on

doCompilerCheck = yes

ifeq ($(CHECK),off)
	doCompilerCheck = no
endif

ifeq ($(CHECK),OFF)
	doCompilerCheck = no
endif


# work out what compiler we're using
# options are:
#
#    gcc
#    g++
#    clang
#

# we only check major and minor version numbers to
# determine if the version of the compiler is good

# set major and minor version numbers required

# Gcc, g++ 7.2.x
gnuVMajorReq = 7
gnuVMinorReq = 2

# clang 7.0.x
clangVMajorReq = 7
clangVMinorReq = 0

CPPVersion = unknown

isGNU = no
isCLANG = no
isCompilerOK = no

ifeq ($(doCompilerCheck),yes)
	# compiler checks are on

	# check compiler

	ifeq ($(CPP),g++)
		isGNU = yes
	else
		ifeq ($(CPP),gcc)
			isGNU = yes
		else
			ifeq ($(CPP),clang)
				isCLANG = yes
			endif
		endif
	endif


	ifeq ($(isGNU),yes)
		# check compiler version
		CPPVersion := $(shell $(CPP) --version | sed -n "s/^.*(GCC)\s*\([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p")
		gnuVMajor_EQ := $(shell expr `echo $(CPPVersion) | cut -f1 -d.` \== $(gnuVMajorReq))
		gnuVMinor_GTEQ := $(shell expr `echo $(CPPVersion) | cut -f2 -d.` \>= $(gnuVMinorReq))
		gnuVMajor_GT := $(shell expr `echo $(CPPVersion) | cut -f1 -d.` \> $(gnuVMajorReq))

		ifeq ($(gnuVMajor_GT),1)
			isCompilerOK = yes
		else
			ifeq ($(gnuVMajor_EQ),1)
				ifeq ($(gnuVMinor_GTEQ), 1)
					isCompilerOK = yes
				endif
			endif
		endif
	else
		ifeq ($(isCLANG),yes)
			# check compiler version
			CPPVersion := $(shell $(CPP) --version | sed -n "s/^.*version\s*\([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p")
			clangVMajor_EQ := $(shell expr `echo $(CPPVersion) | cut -f1 -d.` \== $(clangVMajorReq))
			clangVMinor_GTEQ := $(shell expr `echo $(CPPVersion) | cut -f2 -d.` \>= $(clangVMinorReq))
			clangVMajor_GT := $(shell expr `echo $(CPPVersion) | cut -f1 -d.` \> $(clangVMajorReq))

			ifeq ($(clangVMajor_GT),1)
				isCompilerOK = yes
			else
				ifeq ($(clangVMajor_EQ),1)
					ifeq ($(clangVMinor_GTEQ), 1)
						isCompilerOK = yes
					endif
				endif
			endif
		endif
	endif
else
	# compiler checkes are off
	isCompilerOK = yes
endif

# if compiler version is good, compile, otherwise complain
ifeq ($(isCompilerOK),no)
BUILD:
# bad compiler or version - announce error

# if make wasn't so finicky I could indent these...
ifeq ($(isGNU),yes)
	$(info GCC/g++ version $(gnuVMajorReq).$(gnuVMinorReq).0 or greater is required to compile COMPAS)
else

ifeq ($(isCLANG),yes)
	$(info clang version $(clangVMajorReq).$(clangVMinorReq).0 or greater is required to compile COMPAS)
else
	$(info COMPAS must be compiled with GNU gcc, GNU g++ or clang)
endif

endif

else
# compiler is good - build COMPAS

$(info Building COMPAS with $(CPP) version $(CPPVersion))

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


endif

