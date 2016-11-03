###### User-configurable parameteres

# CPLEX parameters
cplex_home=$(HOME)/sys/opt/optim/ibm/ILOG/CPLEX_Studio1263

# dcsxx-commons parameters
dcsxx_commons_home=$(HOME)/Projects/src/dcsxx-commons

# gtpack parameters
gtpack_home=$(HOME)/Projects/src/gtpack

# Boost parameters
boost_home=$(HOME)/sys/src/git/boost
################################################################################


# Variables
cplex_incs=-I$(cplex_home)/cplex/include -I$(cplex_home)/concert/include -DIL_STD
cplex_libflags=-L$(cplex_home)/cplex/lib/x86-64_linux/static_pic -L$(cplex_home)/concert/lib/x86-64_linux/static_pic
cplex_libs=-lilocplex -lconcert -lcplex -lm -lpthread


# general settings
CXXFLAGS+=-Wall -Wextra -ansi -pedantic
CXXFLAGS+=-g

# dcsxx-commons
CXXFLAGS+=-I$(dcsxx_commons_home)/inc

# gtpack
CXXFLAGS+=-I$(gtpack_home)/include

# Boost
CXXFLAGS+=-I$(boost_home)

# dcsxx-cloud-gt
CXXFLAGS+=-I./inc
CXXFLAGS+=-DDCS_CLOUD_GT_HAVE_CPLEX_SOLVER
CXXFLAGS+=$(cplex_incs)
LDFLAGS+=$(cplex_libflags)
LDLIBS+=$(cplex_libs)
LDLIBS+=-lm


.PHONY: all clean


all: src/sim

clean:
	$(RM) src/sim.o src/sim
