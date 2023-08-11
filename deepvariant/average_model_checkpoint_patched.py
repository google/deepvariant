# Copyright 2023 Google LLC.
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
"""Override _save_model in AverageModelCheckpoint to fix a bug (internal)."""

# TODO: Externally, we can try to use older version of
# average_model_checkpoint.py. For example: r0.13, so that it will be compatible
# to this version we're patching internally.

from tensorflow_addons.callbacks.average_model_checkpoint import AverageModelCheckpoint
from tensorflow_addons.optimizers.average_wrapper import AveragedOptimizerWrapper


# pylint: disable=bad-super-call
class AverageModelCheckpointPatched(AverageModelCheckpoint):
  """_save_model from AverageModelCheckpoint is missing the batch arg."""

  def __init__(
      self,
      update_weights: bool,
      filepath: str,
      monitor: str = 'val_loss',
      verbose: int = 0,
      save_best_only: bool = False,
      save_weights_only: bool = False,
      mode: str = 'auto',
      save_freq: str = 'epoch',
      best_metrics_from_checkpoint: float = float('-inf'),
      **kwargs,
  ):
    super().__init__(
        update_weights,
        filepath,
        monitor,
        verbose,
        save_best_only,
        save_weights_only,
        mode,
        save_freq,
        **kwargs,
    )
    if (
        best_metrics_from_checkpoint is not None
        and best_metrics_from_checkpoint > float('-inf')
    ):
      self.best = best_metrics_from_checkpoint

  def _save_model(self, epoch, batch, logs):
    optimizer = self._get_optimizer()
    assert isinstance(optimizer, AveragedOptimizerWrapper)

    if self.update_weights:
      optimizer.assign_average_vars(self.model.variables)

      # Use super call two levels up to bypass the problematic function
      return super(AverageModelCheckpoint, self)._save_model(epoch, batch, logs)
    else:
      # Note: `model.get_weights()` gives us the weights (non-ref)
      # whereas `model.variables` returns references to the variables.
      non_avg_weights = self.model.get_weights()
      optimizer.assign_average_vars(self.model.variables)
      # result is currently None, since `super._save_model` doesn't
      # return anything, but this may change in the future.

      # Use call two levels up to bypass the problematic function
      result = super(AverageModelCheckpoint, self)._save_model(
          epoch, batch, logs
      )
      self.model.set_weights(non_avg_weights)
      return result

  def on_test_begin(self, logs):
    if self.update_weights:
      self._get_optimizer().swap_weights()

  def on_test_end(self, logs):
    if self.update_weights:
      self._get_optimizer().swap_weights()
