# Copyright 2017 Google Inc.
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
by variantutils and variantcallutils rather than being used directly.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import types

from deepvariant.util.genomics import struct_pb2


def _add_field_with_type(field_map, field_name, value, value_type):
  """Adds values to a particular map field containing a ListValue."""
  if not isinstance(value, (list, types.GeneratorType, tuple)):
    value = [value]
  struct_values = [
      struct_pb2.Value(**{value_type + '_value': v}) for v in value
  ]
  field_map[field_name].values.extend(struct_values)


def _set_field_with_type(field_map, field_name, value, value_type):
  """Sets values to a particular map field containing a ListValue."""
  if field_name in field_map:
    del field_map[field_name]
  _add_field_with_type(field_map, field_name, value, value_type)


def add_number_field(field_map, field_name, value):
  """Appends the given number value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The number value(s) to append to the field. This can be a single
      number or a list of numbers.
  """
  _add_field_with_type(field_map, field_name, value, 'number')


def set_number_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given number value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The number value(s) to set the field to. This can be a single number
      or a list of numbers.
  """
  _set_field_with_type(field_map, field_name, value, 'number')


def add_int_field(field_map, field_name, value):
  """Appends the given int value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The int value(s) to append to the field. This can be a single
      int or a list of ints.
  """
  _add_field_with_type(field_map, field_name, value, 'int')


def set_int_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given int value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The int value(s) to set the field to. This can be a single int
      or a list of ints.
  """
  _set_field_with_type(field_map, field_name, value, 'int')


def add_string_field(field_map, field_name, value):
  """Appends the given string value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The string value(s) to append to the field. This can be a single
      string or a list of strings.
  """
  _add_field_with_type(field_map, field_name, value, 'string')


def set_string_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given string value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The int value(s) to set the field to. This can be a single string or
      a list of strings.
  """
  _set_field_with_type(field_map, field_name, value, 'string')


def add_bool_field(field_map, field_name, value):
  """Appends the given boolean value(s) to field_map[field_name].

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to append value to.
    value: The boolean value(s) to append to the field. This can be a single
      boolean or a list of booleans.
  """
  _add_field_with_type(field_map, field_name, value, 'bool')


def set_bool_field(field_map, field_name, value):
  """Sets field_map[field_name] with the given boolean value(s).

  Args:
    field_map: Map(str --> ListValue) to modify.
    field_name: str. The name of the field to set.
    value: The boolean value(s) to set the field to. This can be a single
      boolean or a list of booleans.
  """
  _set_field_with_type(field_map, field_name, value, 'bool')
