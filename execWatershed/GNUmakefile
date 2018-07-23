# -*- Mode: Makefile -*-

DIM=3

# trace the chain of include makefiles
makefiles += example_watershed

# This is an embedded boundary example
#USE_PETSC=TRUE
#petscincflags = -I/$(PETSC_DIR)/include -I/$(PETSC_DIR)/$(PETSC_ARCH)/include
#petsclibflags = -L/$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

USE_EB=TRUE

# the base name(s) of the application(s) in this directory
ebase = watershedDriver

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../../Chombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic  AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../srcSubsurface ../srcWatershed ../muparser ../../NewPetscSolver

# input file for 'run' target
INPUT = example.hillslope.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
syslibflags +=-L/export/users/mookie/valgrind_parallel/valgrind-3.5.0/mpi/ -lmpiwrap-amd64-linux -Xlinker
