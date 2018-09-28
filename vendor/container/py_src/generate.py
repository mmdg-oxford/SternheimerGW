# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import textwrap

def fortran(file_description):
  """ Generate Fortran file for reading/writing container.

  file_description -- A dictionary that contains a 'name' variable specifying
    the generated module as well as a set of variables the container should
    contain. Each of them is expected to contain a 'type' and a 'dim' variable
    that define the type and the dimensionality of the generated arrays.
  """
  container = license_info()
  container += module_header(file_description)
  container += module_body(file_description)
  container += module_footer(file_description['name'])
  return container

def license_info():
  return textwrap.dedent("""
    ! This Source Code Form is subject to the terms of the Mozilla Public
    ! License, v. 2.0. If a copy of the MPL was not distributed with this
    ! file, You can obtain one at http://mozilla.org/MPL/2.0/.
  """.strip('\n'))

def module_header(file_description):
  type_name = '{0}_container'.format(file_description['name'])
  type_string = ''
  for key in file_description:
      if key == 'name': continue
      dim_string = "(:" + ",:" * (file_description[key]['dim'] - 1) + ")"
      type_string += "\n    {0}, ALLOCATABLE :: {1}".format(
        file_description[key]['type'].upper(), key + dim_string)
  var_name = ['var_{0}'.format(key) for key in file_description if key != 'name']
  public_string = ', '.join(var_name + [type_name])
  enum_string = 'ENUMERATOR :: {0} = 1'.format(var_name[0])
  if len(var_name) > 1:
      enum_string += '\n    ENUMERATOR ' + ', '.join(var_name[1:])
  return """MODULE {0}

  USE container_interface

  IMPLICIT NONE

  PRIVATE
  PUBLIC {1}

  ENUM, BIND(C)
    {2}
  END ENUM

  TYPE, EXTENDS(container_type) :: {3}{4}
  CONTAINS
    PROCEDURE :: internal_init => init
    PROCEDURE :: internal_read_variable => read_variable
    PROCEDURE :: internal_read_element => read_element
    PROCEDURE :: internal_write_variable => write_variable
    PROCEDURE :: internal_write_element => write_element
    PROCEDURE :: internal_update_offset => update_offset
  END TYPE {3}

""".format(file_description['name'], public_string, enum_string, type_name, type_string)

  type_string = """  TYPE {0}
    INTEGER filehandle
    LOGICAL :: valid = .FALSE.
""".format(type_name)
  for key in file_description:
    if key == 'name': continue
    type_string += "    " + file_description[key]['type'].upper() + ", ALLOCATABLE :: " + key + "(:" + ",:" * (file_description[key]['dim'] - 1) + ")\n"
  type_string += "  END TYPE " + type_name + "\n\n"
  return type_string

def module_body(file_description):
  body = "CONTAINS\n"
  body += init(file_description)
  body += write_variable(file_description)
  body += write_element(file_description)
  body += read_variable(file_description)
  body += read_element(file_description)
  body += update_offset(file_description)
  return body

def module_footer(name):
  return "\nEND MODULE {0}\n".format(name)

def init(file_description):
  dim = [file_description[key]['dim'] for key in file_description if key != 'name']
  return """
  SUBROUTINE init(container)
    !
    CLASS({0}_container), INTENT(OUT) :: container
    container%num_dim = {1}
    !
  END SUBROUTINE init
""".format(file_description['name'], dim)

def write_variable(file_description):
  write_string = ''
  for key in file_description:
    if key == 'name': continue
    dim = file_description[key]['dim']
    mpi_type = _mpi_type(file_description[key]['type'])
    write_string += """
    CASE (var_{0})
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%{0}), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%{0}), {1}, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%{0}, SIZE(container%{0}), &
        {2}, MPI_STATUS_IGNORE, ierr)""".format(key, dim, mpi_type)
  return """
  SUBROUTINE write_variable(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS({0}_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    !
    SELECT CASE (config%variable){1}
    END SELECT
    !
  END SUBROUTINE write_variable
""".format(file_description['name'], write_string)

def write_element(file_description):
  write_string = ''
  for key in file_description:
    if key == 'name': continue
    dim = file_description[key]['dim']
    mpi_type = _mpi_type(file_description[key]['type'])
    write_string += """
    CASE (var_{0})
      CALL container%increase_offset(dims, {1}, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%{0}, SIZE(container%{0}), &
        {1}, MPI_STATUS_IGNORE, ierr)""".format(key, mpi_type)
  return """
  SUBROUTINE write_element(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS({0}_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable){1}
    END SELECT
    !
  END SUBROUTINE write_element
""".format(file_description['name'], write_string)

def read_variable(file_description):
  read_string = ''
  for key in file_description:
    if key == 'name': continue
    dim_string = _dim_string_variable(file_description[key]['dim'])
    mpi_type = _mpi_type(file_description[key]['type'])
    read_string += """
    CASE (var_{0})
      IF (ALLOCATED(container%{0})) DEALLOCATE(container%{0})
      ALLOCATE(container%{0}({1}))
      CALL mpi_func(container%filehandle, container%{0}, SIZE(container%{0}), &
        {2}, MPI_STATUS_IGNORE, ierr)""".format(key, dim_string, mpi_type)
  return """
  SUBROUTINE read_variable(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS({0}_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    SELECT CASE (config%variable){1}
    END SELECT
    !
  END SUBROUTINE read_variable
""".format(file_description['name'], read_string)

def _dim_string_variable(dim):
  dim_string = ''
  for i in range(dim):
    dim_string += 'dims({0}), '.format(i + 1)
  return dim_string.strip(', ')

def read_element(file_description):
  read_string = ''
  for key in file_description:
    if key == 'name': continue
    dim_string = _dim_string_element(file_description[key]['dim'])
    mpi_type = _mpi_type(file_description[key]['type'])
    read_string += """
    CASE (var_{0})
      CALL container%increase_offset(dims, {2}, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%{0})) DEALLOCATE(container%{0})
      ALLOCATE(container%{0}({1}))
      CALL mpi_func(container%filehandle, container%{0}, SIZE(container%{0}), &
        {2}, MPI_STATUS_IGNORE, ierr)""".format(key, dim_string, mpi_type)
  return """
  SUBROUTINE read_element(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS({0}_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable){1}
    END SELECT
    !
  END SUBROUTINE read_element
""".format(file_description['name'], read_string)

def update_offset(file_description):
  offset_string = ''
  for key in file_description:
    if key == 'name': continue
    mpi_type = _mpi_type(file_description[key]['type'])
    offset_string += """
    CASE (var_{0})
      CALL container%increase_offset(dimension, {1}, ierr)""".format(key, mpi_type)
  return """
  SUBROUTINE update_offset(container, variable, dimension, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS({0}_container), INTENT(INOUT) :: container
    INTEGER, INTENT(IN) :: variable, dimension(:)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER(KIND=MPI_OFFSET_KIND) offset
    !
    SELECT CASE (variable){1}
    END SELECT
    IF (ierr /= no_error) RETURN
    CALL container%update_offset(ierr)
    !
  END SUBROUTINE update_offset
""".format(file_description['name'], offset_string)

def _dim_string_element(dim):
  dim_string = ''
  for i in range(dim - 1):
    dim_string += 'dims({0}), '.format(i + 1)
  return dim_string + "1"

def _mpi_type(fortran_type):
  type_dict = {
    'integer': 'MPI_INTEGER',
    'real': 'MPI_REAL',
    'real(dp)': 'MPI_DOUBLE_PRECISION',
    'complex': 'MPI_COMPLEX',
    'complex(dp)': 'MPI_DOUBLE_COMPLEX',
  }
  return type_dict[fortran_type.lower()]
