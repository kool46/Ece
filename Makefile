# Executable name
EXE=spece

# Objects list
OBJ=kinds.o constants.o specelib.o cubic.o intervals.o phys_constants.o spece.o spece_data.o utils.o

# Alternative search paths
vpath %.f90 src
vpath %.mod bin
vpath $(EXE) bin
vpath %.o bin

##
## Git information (used in the version string)
##

# Short hash of the latest commit
GIT_REV   ?= $(shell git rev-parse --short HEAD)

# Whether the worktree and the latest commit differs
GIT_DIRTY ?= $(shell test -n "$$(git status --porcelain)" && echo "-dirty")

##
## Fortran compiler and flags
##

# Fortran compiler name and flags
FC=gfortran
FFLAGS=-O3 -cpp -DREVISION=\"$(GIT_REV)$(GIT_DIRTY)\"

all: $(EXE)

# Build executable from object files
$(EXE): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# Dependencies on modules
spece.o: constants.mod cubic.mod intervals.mod kinds.mod phys_constants.mod spece_data.mod utils.mod
spece_data.o: kinds.mod
specelib.o: constants.mod kinds.mod
cubic.o: constants.mod kinds.mod
intervals.o: constants.mod kinds.mod
phys_constants.o: constants.mod kinds.mod
constants.o: kinds.mod
utils.o: kinds.mod

# Check existence of module information files
%.mod: %.o
	@if [ ! -f $@ ]; then \
	  rm $<; \
	  $(MAKE) $<; \
	fi

# General object compilation command
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

install:
	@if [ -f $(EXE) ]; then \
	  mv $(EXE) $(HOME)/bin/; \
	else \
	  echo File $(EXE) does not exist. Run \'make\' first; \
	fi

.PHONY: clean
# Remove output files
clean:
	rm -rf *.o *.mod $(EXE)
