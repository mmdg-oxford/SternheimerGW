# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
LIB = libcontainer.a
MPIF90 = mpif90
AR = ar
ARFLAGS = rv
RANLIB = ranlib

lib: $(LIB)

$(LIB) : $(OBJ)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean:
	rm -f $(OBJ) $(LIB)

%.o: %.f90
	$(MPIF90) -c $(FFLAGS) $<

-include make.inc
-include make.depend
