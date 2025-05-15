#!/bin/bash
# Copyright 2025 Google LLC.
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

# Runs the given command in a dry run mode, printing the command to be executed
# to stdout and the command's exit status to stderr.
#
# Usage:
#   run command [--skip-dry-run-print]
#
# Options:
#   --skip-dry-run-print: Disables the output in the dry run mode.
#
function run() {
  local no_dry_echo=false
  local actual_command_args=("$@") # Copy arguments to a new array

  # Check for our special flag as the first argument
  if [[ "${1:-}" == "--skip-dry-run-print" ]]; then
    no_dry_echo=true
    # Remove the flag from the arguments to be executed/printed
    actual_command_args=("${@:2}")
  fi

  if [[ "${DRY_RUN:-false}" == "true" ]]; then
    if [[ "$no_dry_echo" == "false" ]]; then
      # Prints out command to stdout and a [DRY RUN] tag to stderr.
      1>&2 printf "[DRY RUN] " && echo "${actual_command_args[*]}"
    # else: do nothing, effectively suppressing the echo for this specific command
    fi
    return 0 # Dry run itself is considered a success.
  else
    # Print the command to be executed (mimics 'set -x' a bit)
    # We use the original "$*" here if we want to see --skip-dry-run-print if it was passed
    # or actual_command_args if we don't want to see it in the trace.
    # Let's print the actual command that will be run:
    echo "+ ${actual_command_args[*]}"
    eval "${actual_command_args[@]}" # Use array expansion for safety with eval
    local exit_status=$?
    return "$exit_status"
  fi
}
