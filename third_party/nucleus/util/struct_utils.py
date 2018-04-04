# Copyright 2018 Google Inc.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""Struct proto utilities.

This class provides wrappers for conveniently interacting with protos defined
in struct.proto, mostly ListValue and Value objects. It should primarily be used
by variant_utils and variantcallutils rather than being used directly.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import types

from third_party.nucleus.protos import struct_pb2

# Field names of values defined in struct_pb2.Value.
_BOOL_TYPE = 'bool_value'
_INT_TYPE = 'int_value'
_NUMBER_TYPE = 'number_value'
_STRING_TYPE = 'string_value'


def _add_field_with_type(field_map, field_name, value, value_type):
  """Adds values to a particular map field containing a ListValue."""
  if not isinstance(value, (list, types.GeneratorType, tuple)):
    value = [value]
  struct_values = [struct_pb2.Value(**{value_type: v}) for v in value]
  field_map[field_name].values.extend(struct_values)


def _set_field_with_type(field_map, field_name, value, value_type):
  """Sets values to a particular map field containing a ListValue."""
  if field_name in field_map:
    del field_map[field_name]
  _add_field_with_type(field_map, field_name, value, value_type)


def _get_field_with_type(field_map, field_name, is_single_field, value_type):
  fields = [getattr(v, value_type) for v in field_map[field_name].values]
  return fields[0] if is_single_field and fields else fields


def add_number_field(field_map, field_name, value):
  """Appends the given number value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The number value(s) to append to the field. This can be a single
      number or a list of numbers.
  """
  _add_field_with_type(field_map, field_name, value, _NUMBER_TYPE)


def set_number_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given number value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The number value(s) to set the field to. This can be a single number
      or a list of numbers.
  """
  _set_field_with_type(field_map, field_name, value, _NUMBER_TYPE)


def get_number_field(field_map, field_name, is_single_field=False):
  """Returns the number value(s) stored in `field_map[field_name]`.

  If the field_name is not present in field_map, the empty list is returned.

  Args:
    field_map: Map(str --> ListValue) of interest.
    field_name: str. The name of the field to extract number values from.
    is_single_field: bool. If True, return the first number value stored (it
      should be the only one in the field). Otherwise, return the list of
      numbers.

  Returns:
    The number value(s) stored in the field_map under this field_name.
  """
  return _get_field_with_type(field_map, field_name, is_single_field,
                              _NUMBER_TYPE)


def add_int_field(field_map, field_name, value):
  """Appends the given int value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The int value(s) to append to the field. This can be a single
      int or a list of ints.
  """
  _add_field_with_type(field_map, field_name, value, _INT_TYPE)


def set_int_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given int value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The int value(s) to set the field to. This can be a single int
      or a list of ints.
  """
  _set_field_with_type(field_map, field_name, value, _INT_TYPE)


def get_int_field(field_map, field_name, is_single_field=False):
  """Returns the int value(s) stored in `field_map[field_name]`.

  If the field_name is not present in field_map, the empty list is returned.

  Args:
    field_map: Map(str --> ListValue) of interest.
    field_name: str. The name of the field to extract int values from.
    is_single_field: bool. If True, return the first int value stored (it
      should be the only one in the field). Otherwise, return the list of
      ints.

  Returns:
    The int value(s) stored in the field_map under this field_name.
  """
  return _get_field_with_type(field_map, field_name, is_single_field, _INT_TYPE)


def add_string_field(field_map, field_name, value):
  """Appends the given string value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The string value(s) to append to the field. This can be a single
      string or a list of strings.
  """
  _add_field_with_type(field_map, field_name, value, _STRING_TYPE)


def set_string_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given string value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The int value(s) to set the field to. This can be a single string or
      a list of strings.
  """
  _set_field_with_type(field_map, field_name, value, _STRING_TYPE)


def get_string_field(field_map, field_name, is_single_field=False):
  """Returns the string value(s) stored in `field_map[field_name]`.

  If the field_name is not present in field_map, the empty list is returned.

  Args:
    field_map: Map(str --> ListValue) of interest.
    field_name: str. The name of the field to extract string values from.
    is_single_field: bool. If True, return the first string value stored (it
      should be the only one in the field). Otherwise, return the list of
      strings.

  Returns:
    The string value(s) stored in the field_map under this field_name.
  """
  return _get_field_with_type(field_map, field_name, is_single_field,
                              _STRING_TYPE)


def add_bool_field(field_map, field_name, value):
  """Appends the given boolean value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The boolean value(s) to append to the field. This can be a single
      boolean or a list of booleans.
  """
  _add_field_with_type(field_map, field_name, value, _BOOL_TYPE)


def set_bool_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given boolean value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The boolean value(s) to set the field to. This can be a single
      boolean or a list of booleans.
  """
  _set_field_with_type(field_map, field_name, value, _BOOL_TYPE)


def get_bool_field(field_map, field_name, is_single_field=False):
  """Returns the bool value(s) stored in `field_map[field_name]`.

  If the field_name is not present in field_map, the empty list is returned.

  Args:
    field_map: Map(str --> ListValue) of interest.
    field_name: str. The name of the field to extract bool values from.
    is_single_field: bool. If True, return the first bool value stored (it
      should be the only one in the field). Otherwise, return the list of
      bools.

  Returns:
    The bool value(s) stored in the field_map under this field_name.
  """
  return _get_field_with_type(field_map, field_name, is_single_field,
                              _BOOL_TYPE)
