import sys
sys.path.append('../../../vendor/container/py_src')
import generate

file_description = {
  'name': 'gw_data',
  'coul': {
    'type': 'complex(dp)',
    'dim': 4,
  },
  'corr': {
    'type': 'complex(dp)',
    'dim': 4,
  },
  'exch': {
    'type': 'complex(dp)',
    'dim': 3,
  },
  'k_point': {
    'type': 'real(dp)',
    'dim': 2,
  },
}

container_module = generate.fortran(file_description)
print container_module
