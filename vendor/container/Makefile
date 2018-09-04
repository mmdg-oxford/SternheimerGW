# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
lib:
	make -C f90_src

test:
	make -C test
	make -C py_src

clean:
	make -C f90_src clean
	make -C example clean
	make -C test clean

.PHONY: lib test clean
