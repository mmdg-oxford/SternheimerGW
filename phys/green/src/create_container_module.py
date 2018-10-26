import sys
sys.path.append('../../../vendor/container/py_src')
import generate

file_description = {
  'name': 'linear_problem',
  'omega': {
    'type': 'complex(dp)',
    'dim': 1,
  },
  'hamil': {
    'type': 'complex(dp)',
    'dim': 2,
  },
  'rhs': {
    'type': 'complex(dp)',
    'dim': 1,
  },
  'green': {
    'type': 'complex(dp)',
    'dim': 2,
  },
}

container_module = generate.fortran(file_description)
print container_module
