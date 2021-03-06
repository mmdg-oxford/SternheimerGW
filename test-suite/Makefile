#------------------------------------------------------------------------------
#
# This file is part of the SternheimerGW code.
# Parts of this file are taken from the Quantum ESPRESSO software
# P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
#
# Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
# Henry Lambert, Martin Schlipf, and Feliciano Giustino
#
# SternheimerGW is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SternheimerGW is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SternheimerGW. If not, see
# http://www.gnu.org/licenses/gpl.html .
#
#------------------------------------------------------------------------------

include ./ENVIRONMENT

default :
	@echo 'To run Sternheimer GW test-suite, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' run-tests          : run-tests-serial (default)'
	@echo ' run-tests-serial   : run all tests in SERIAL' 
	@echo ' run-tests-parallel : run all tests in PARALLEL (4 MPI)'
	@echo ' compare            : compare last output generated with reference'
	@echo ' clean              : clean output (our/err files) of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

pseudo :
	@./check_pseudo.sh gw_

run-tests : run-tests-serial

run-tests-serial : run-tests-gw-serial

run-tests-parallel : run-tests-gw-parallel


run-tests-gw-serial : pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=gw_all

run-tests-gw-parallel : pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=gw_all

run-custom-test : run-custom-test-serial

run-custom-test-serial : pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-custom-test-parallel : pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi


create-reference-gw :
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=gw_all make-benchmarks

create-custom-reference :
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=$(testdir) make-benchmarks

compare : compare-gw

compare-gw :
	${TESTCODE_DIR}/bin/testcode.py --category=gw_all  --verbose compare 2>&1 | tee out.SternheimerGW.`date +%Y%m%d_%H%M%S`

compare-custom-gw :
	${TESTCODE_DIR}/bin/testcode.py --category=$(testdir) --verbose compare 2>&1 | tee out.$(testdir).`date +%Y%m%d_%H%M%S`

clean:
	@for x in `find gw_* -name "test*"`; do rm -rf $$x; done
	@for x in `ls -d gw_*`; do rm -rf $$x/tmp; rm -f $$x/fort.1000; done
	@rm -rf out.SternheimerGW.*

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find gw_* -name "benchmark*"`; do rm -rf $$x; done
