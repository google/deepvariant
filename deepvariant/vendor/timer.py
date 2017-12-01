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
"""This module implements classes useful for timing how long an action took.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import time


class TimerError(Exception):
  """Timer Error Exception."""


class Timer(object):
  """Measure time elapsed since some event.

  Example usage:

  from deepvariant.vendor import timer
  a_timer = timer.Timer()
  a_timer.Start()

  # Do stuff here ...

  a_timer.Stop()
  if a_timer.GetDuration() > 100:
    print 'This took too long'

  Another way to use this class is the with statement:

  with timer.Timer() as t:
    # Do time consuming stuff ...
  print 'Time consuming stuff took %f seconds.' % t.GetDuration()
  """

  def __init__(self):
    """Initializes a timer."""
    self.duration = 0
    self.running = False
    self.start = None

  def Start(self):
    """Resets and starts a timer."""

    self._StartInternal()

  def Stop(self):
    """Stops a timer, and records the duration. Stop is idempotent.

    Returns:
      The duration, i.e., the time (in seconds) for which the timer ran.
    """
    self._StopInternal()
    return self.duration

  def __enter__(self):
    """Resets and starts a timer.

    This allows Timer to be used as a ContextManager type.

    Returns:
      The object itself so that it can be bound to a variable in a with
      statement.
    """
    self._StartInternal()
    return self

  def __exit__(self, unused_ex_type, unused_ex, unused_ex_trace):
    """Stops a timer and records the duration.

    This allows Timer to be used as a ContextManager type.

    Returns:
      False.  This means that any exceptions raised within the body of the
      caller's with statement will propagate up the stack without interference.
    """
    self._StopInternal()
    return False

  def _StartInternal(self):
    """Resets and starts a timer.

    Common implementation for Start and __enter__.
    """
    self.start = time.time()
    self.running = 1

  def _StopInternal(self):
    """Stops a timer and records the duration.

    Common implementation for Stop and __exit__.
    """
    if self.running:
      self.duration = time.time() - self.start
      self.running = 0

  def IsRunning(self):
    return self.running

  def __str__(self):
    return str(self.GetDuration())

  def GetDuration(self):
    """Returns the elapsed time, in seconds.

    Returns:
       If timer is still running: the time (in seconds) elapsed since the start
       Otherwise: the time (in seconds) for which the timer ran.
    """
    if self.running:
      return time.time() - self.start
    else:
      return self.duration

  def GetStartTime(self):
    """Returns start time of the timer.

    Returns:
       The start time of the timer, floating-point seconds since the epoch.
    Raises:
       TimerError: if the timer was not started
    """
    if self.start:
      return self.start
    else:
      raise TimerError('TimerNotStarted')

  def GetStopTime(self):
    """Returns stop time of the timer.

    Returns:
       The stop time of the timer, floating-point seconds since the epoch.
    Raises:
       TimerError: if the timer was never started or is still running
    """
    if not self.start:
      raise TimerError('TimerNotStarted')
    elif self.running:
      raise TimerError('TimerStillRunning')
    else:
      return self.start + self.duration


class TimerStart(Timer):
  """A timer that automatically starts when you construct it.

  This is otherwise identical in interface to the Timer class in this module.
  """

  def __init__(self):
    Timer.__init__(self)
    self.Start()


class MultiIntervalTimer(Timer):
  """A timer that records cumulative time intervals.

  Example usage:

  from deepvariant.vendor import timer
  a_timer = timer.MultiIntervalTimer()

  # Do stuff the first time
  a_timer.Start()
  # ... stuff ...
  a_timer.Stop()

  # Do stuff the second time (repeat as many times as you like)
  a_timer.Start()
  # ... stuff ...
  a_timer.Stop()

  if a_timer.GetDuration() > 1000:
    print 'Total time spent doing stuff is too much!'

  Another way to use this class is the with statement:

  t = timer.MultiIntervalTimer()
  with t:
    # Do stuff the first time

  with t:
    # Do stuff the second time (repeat as many times as you like)

  print 'Total time spent doing stuff was %f seconds.' % t.GetDuration()
  """

  def __init__(self):
    super(MultiIntervalTimer, self).__init__()
    self.stop = None

  def Reset(self):
    """Resets the measured cumulative time to zero."""
    self.duration = 0

  def _StopInternal(self):
    """Stops a timer and adds the current duration to the total.

    Common implementation for Stop and __exit__.
    """
    if self.running:
      self.stop = time.time()
      self.duration += self.stop - self.start
      self.running = 0

  def GetDuration(self):
    """Returns the total elapsed time, in seconds.

    Returns:
       The total time (in seconds) accumulated so far, regardless of counter
       state (i.e. if stopped, total time spent in all running windows; if
       running, total time spent in all previous running windows plus time
       elapsed from the last call to Start() until GetDuration() was called).
    """
    if self.running:
      return self.duration + (time.time() - self.start)
    else:
      return self.duration

  def GetStopTime(self):
    """Returns the last stop time of the timer.

    Returns:
       The last stop time of the timer, same format as time.time()
    Raises:
       TimerError: if the timer was never started or is still running
    """
    if not self.start:
      return super(MultiIntervalTimer, self).GetStopTime()
    elif self.running:
      return super(MultiIntervalTimer, self).GetStopTime()
    else:
      return self.stop
