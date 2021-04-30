# LIFE makefile

# Compiler command
CC=g++
CFLAGS=-O3 -std=c++0x -fopenmp -Wall -Wextra

# Executable
EXE=LIFE

# Location of source, header and object files
DIR=.
SDIR=$(DIR)/src
HDIR=$(DIR)/inc
ODIR=$(DIR)/obj

# Get the sources and object files
LIFESRCS:=FEMBody.cpp FEMElement.cpp FEMNode.cpp Grid.cpp GridUtils.cpp IBMBody.cpp IBMNode.cpp IBMSupport.cpp Objects.cpp Utils.cpp main.cpp
LIFEOBJS:=$(addprefix $(ODIR)/,$(notdir $(LIFESRCS:.cpp=.o)))
LIFEHDRS:=$(wildcard $(HDIR)/*.h)

TESTFEMSRCS:=TestFEM.cpp
TESTFEMOBJS:=$(addprefix $(ODIR)/,$(notdir $(TESTFEMSRCS:.cpp=.o)))
FEMLIBOBJS:=$(addprefix $(ODIR)/, FEMBody.o FEMElement.o FEMNode.o Utils.o)

ifeq ($(LAPACK_LIB),)
LAPACK_LIB := -llapack
endif

OBJS:=$(LIFEOBJS) $(TESTFEMOBJS)
# Include and library files
INC=
LIB=$(LAPACK_LIB) -lboost_system -lboost_filesystem

-include make.config

# Build LIFE
$(EXE): $(LIFEOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIB)

# Build object files
$(OBJS): $(ODIR)/%.o : $(SDIR)/%.cpp $(LIFEHDRS)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Clean the project
.PHONY: clean
clean:
	rm -rf $(EXE) TestFEM $(ODIR) Results *.out && mkdir $(ODIR)

# FEM test program
TestFEM: $(TESTFEMOBJS) $(ODIR)/libLIFEFEM.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIB)

# Build the FEM library
$(ODIR)/libLIFEFEM.a: $(FEMLIBOBJS)
	$(AR) rcs $@ $^
