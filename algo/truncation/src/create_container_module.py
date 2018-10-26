import sys
sys.path.append('../../../vendor/container/py_src')
import generate

file_description = {
  'name': 'trunc_data',
  'cell': {
    'type': 'real(dp)',
    'dim': 3,
  },
  'cutoff': {
    'type': 'real(dp)',
    'dim': 1,
  },
  'trunc_coul': {
    'type': 'real(dp)',
    'dim': 3,
  },
}

container_module = generate.fortran(file_description)
print container_module
