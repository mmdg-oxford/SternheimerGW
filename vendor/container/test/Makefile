# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

.PHONY: all lib clean
.DEFAULT_GOAL = all

include $(PFUNIT)/include/base.mk

FFLAGS += -I$(PFUNIT)/mod -I../f90_src -I../example

ifeq ($(USEMPI),YES)
   MPIF90 = mpif90
   FPPFLAGS += -DUSE_MPI
endif

EXE = tests$(EXE_EXT)
SRC = $(wildcard *.pf)
OBJ = $(SRC:.pf=$(OBJ_EXT))

LIB = ../f90_src/libcontainer.a ../example/libexample.a

all: lib $(EXE)
	mpirun -np 4 ./$(EXE)

$(EXE): testSuites.inc $(OBJ) $(LIB)
	$(MPIF90) -o $@ -I. $(FFLAGS) $(FPPFLAGS) $(PFUNIT)/include/driver.F90 $(OBJ) -L$(PFUNIT)/lib -lpfunit $(LIB)

%$(OBJ_EXT): %.F90
	$(MPIF90) -c $(FFLAGS) $<

%.F90: %.pf $(LIB)
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

lib:
	make -C ../f90_src lib
	make -C ../example lib

clean:
	rm -f *.io

-include make.depend
