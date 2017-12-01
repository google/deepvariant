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
"""Save and restore flag values after the decorated function completes.

There are many ways to save and restore.  Always use the most convenient method
for a given use case.

Here are examples of each method.  They all call DoStuff() while FLAGS.someflag
is temporarily set to 'foo'.

  # Use a decorator which can override flags via arguments.
  @flagsaver.FlagOverrider(someflag='foo')
  def SomeFunc():
    DoStuff()

  # Use a decorator which does not override flags itself.
  @flagsaver.FlagSaver
  def SomeFunc():
    FLAGS.someflag = 'foo'
    DoStuff()

  # Use a context manager which can optionally override flags via arguments.
  with flagsaver.FlagContext(someflag='foo'):
    DoStuff()

  # Save and restore the flag values yourself.
  saved_flag_values = flagsaver.SaveFlagValues()
  try:
    FLAGS.someflag = 'foo'
    DoStuff()
  finally:
    flagsaver.RestoreFlagValues(saved_flag_values)

We save and restore a shallow copy of each Flag object's __dict__ attribute.
This preserves all attributes of the flag, such as whether or not it was
overridden from its default value.

WARNING: Currently a flag that is saved and then deleted cannot be restored.  An
exception will be raised.  However if you *add* a flag after saving flag values,
and then restore flag values, the added flag will be deleted with no errors.  If
you wish to delete and then restore saved flags, send a CL.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import contextlib
import functools
import inspect


import tensorflow as tf

from tensorflow.python.platform import flags

FLAGS = tf.flags.FLAGS


@contextlib.contextmanager
def FlagContext(**overrides):
  """Context manager that saves flag values on entry and restores them on exit.

  Args:
    **overrides: Set these flag values while within the context manager.

  Yields:
    Nothing, it's a context manager.
  """
  saved_flag_values = SaveFlagValues()
  for name, value in overrides.iteritems():
    setattr(FLAGS, name, value)
  try:
    yield
  finally:
    RestoreFlagValues(saved_flag_values)


def _FlagSaver(func, overrides):
  """Creates a wrapper function that saves/restores flag values.

  Args:
    func: function object - This will be called between saving flags and
        restoring flags.
    overrides: {str: object} - Flag names mapped to their values.  These flags
        will be set after saving the original flag state.

  Returns:
    return value from func()
  """

  @functools.wraps(func)
  def _FlagSaverWrapper(*args, **kwargs):
    """Wrapper function that saves and restores flags."""
    with FlagContext(**overrides):
      return func(*args, **kwargs)

  return _FlagSaverWrapper


def FlagSaver(func):
  """Saves/restores flag values after decorated method completes."""
  if inspect.isclass(func):
    raise TypeError('@flagsaver.FlagSaver cannot be applied to a class.')
  return _FlagSaver(func, {})


class FlagOverrider(object):
  """Overrides flags for the duration of the decorated function call.

  It also restores all original values of flags after decorated method
  completes.
  """

  def __init__(self, **overrides):
    self._overrides = overrides

  def __call__(self, func):
    if inspect.isclass(func):
      raise TypeError('@flagsaver.FlagOverrider cannot be applied to a class.')
    return _FlagSaver(func, self._overrides)


def _CopyFlagDict(flag):
  """Returns a copy of the flag object's __dict__.

  It's mostly a shallow copy of the __dict__, except it also does a shallow
  copy of the validator list.

  Args:
    flag: A flags.Flag instance.

  Returns:
    A copy of the flag object's __dict__.
  """
  copy = flag.__dict__.copy()
  copy['validators'] = list(flag.validators)
  return copy


def SaveFlagValues():
  """Returns copy of flag values as a dict.

  Returns:
    Dictionary mapping keys to values. Keys are flag names, values are
    corresponding __dict__ members. E.g. {'key': value_dict, ...}.
  """
  if hasattr(flags, '_FlagValues'):  # pylint:disable=protected-access
    # In OSS code we use tensorflow/python/platform/flags.py:_FlagValues
    # which is not iterable.
    flag_dict = FLAGS.__dict__['__flags']
    # Make a shallow copy of the flags.
    return {name: flag_dict[name] for name in flag_dict}
  else:
    # FLAGS is iterable and provides __getitem__.
    return {name: _CopyFlagDict(FLAGS[name]) for name in FLAGS}


def RestoreFlagValues(saved_flag_values):
  """Restores flag values based on the dictionary of flag values.

  Args:
    saved_flag_values: {'key': value_dict, ...}
  """
  if hasattr(flags, '_FlagValues'):  # pylint:disable=protected-access
    # In OSS code we use tensorflow/python/platform/flags.py:_FlagValues
    # which is not iterable.
    new_flag_names = FLAGS.__dict__['__flags'].keys()
  else:
    # FLAGS is iterable and provides __getitem__.
    new_flag_names = list(FLAGS)

  for name in new_flag_names:
    value = saved_flag_values.get(name)
    if value is None:
      # If value was not saved delete it
      if hasattr(flags, '_FlagValues'):  # pylint:disable=protected-access
        # In OSS code we use tensorflow/python/platform/flags.py.
        del FLAGS.__dict__['__flags'][name]
      else:
        delattr(FLAGS, name)

    else:
      if hasattr(flags, '_FlagValues'):  # pylint:disable=protected-access
        # In OSS code we use tensorflow/python/platform/flags.py:_FlagValues
        # which is not iterable.
        FLAGS.__dict__['__flags'][name] = value
      else:
        # We bypass the overridden __setitem__ and __setattr__ methods of FLAGS.
        FLAGS[name].__dict__ = value
