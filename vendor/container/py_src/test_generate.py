# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import generate
import pytest
import sys
sys.path.append('../example')
import simple
import two_array
import complex

def test_not_dict():
  with pytest.raises(Exception):
    generate.fortran(0)

def test_simple():
  ref_file = open('../example/simple.f90', 'r')
  ref_content = ref_file.read()
  ref_file.close()
  content = generate.fortran(simple.dict)
  print content
  assert ref_content == content

def test_two_array():
  ref_file = open('../example/two_array.f90', 'r')
  ref_content = ref_file.read()
  ref_file.close()
  content = generate.fortran(two_array.dict)
  print content
  assert ref_content == content

def test_complex_array():
  ref_file = open('../example/complex.f90', 'r')
  ref_content = ref_file.read()
  ref_file.close()
  content = generate.fortran(complex.dict)
  print content
  assert ref_content == content
