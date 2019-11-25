# Copyright 2019 Google LLC.
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

# Lint as: python3
"""Utility functions for visualization and inspection of pileup examples.

Visualization and inspection utility functions enable showing image-like array
data including those used in DeepVariant.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from IPython import display
import numpy as np
from PIL import Image

from third_party.nucleus.protos import variants_pb2


def get_image_array_from_example(example):
  """Decode image/encoded and image/shape of an Example into a numpy array.

  Parse image/encoded and image/shape features from a tensorflow Example and
  decode the image into that shape.

  Args:
    example: a tensorflow Example containing features that include
      "image/encoded" and "image/shape"

  Returns:
    numpy array of dtype np.uint8.
  """
  features = example.features.feature
  img = features["image/encoded"].bytes_list.value[0]
  shape = features["image/shape"].int64_list.value[0:3]
  return np.frombuffer(img, np.uint8).reshape(shape)


def split_3d_array_into_channels(arr):
  """Split 3D array into a list of 2D arrays.

  e.g. given a numpy array of shape (100, 200, 6), return a list of 6 channels,
  each with shape (100, 200).

  Args:
    arr: a 3D numpy array.

  Returns:
    list of 2D numpy arrays.
  """
  return [arr[:, :, i] for i in range(arr.shape[-1])]


def channels_from_example(example):
  """Extract image from an Example and return the list of channels.

  Args:
    example: a tensorflow Example containing features that include
      "image/encoded" and "image/shape"

  Returns:
    list of 2D numpy arrays, one for each channel.
  """
  image = get_image_array_from_example(example)
  return split_3d_array_into_channels(image)


def convert_6_channels_to_rgb(channels):
  """Convert 6-channel image from DeepVariant to RGB for quick visualization.

  The 6 channels are: "read base", "base quality", "mapping quality", "strand",
  "supports variant", "base != reference".

  Args:
    channels: a list of 6 numpy arrays.

  Returns:
    3D numpy array of 3 colors (Red, green, blue).
  """
  base = channels[0]
  # qual is the minimum of base quality and mapping quality at each position
  # 254 is the max value for quality scores because the SAM specification has
  # 255 reserved for unavailable values.
  qual = np.minimum(channels[1], channels[2])
  strand = channels[3]
  # alpha is <supports variant> * <base != reference>
  alpha = np.multiply(channels[4] / 254.0, channels[5] / 254.0)
  return np.multiply(np.stack([base, qual, strand]),
                     alpha).astype(np.uint8).transpose([1, 2, 0])


def adjust_colors_for_png(arr, vmin=0, vmax=255):
  """Scale an array to integers between 0 and 255 to prep it for a PNG image.

  Args:
    arr: numpy array. Input array made up of integers or floats.
    vmin: number. Minimum data value to map to 0. Values below this will be
      clamped to this value and therefore become 0.
    vmax: number. Maximum data value to map to 255. Values above this will be
      clamped to this value and therefore become 255.

  Returns:
    numpy array of dtype np.uint8 (integers between 0 and 255).
  """
  if vmax == 0 or vmax <= vmin:
    raise ValueError("vmin must be non-zero and higher than vmin.")

  # Careful not to modify the original array
  scaled = np.copy(arr)

  # Snap numbers in the array falling outside the range into the range,
  # otherwise they will produce artifacts due to byte overflow
  scaled[scaled > vmax] = vmax
  scaled[scaled < vmin] = vmin

  # Scale the input into the range of vmin to vmax
  if vmin != 0 or vmax != 255:
    scaled = ((scaled - vmin) / (vmax - vmin)) * 255
  return scaled.astype(np.uint8)


def enlarge_image_array(arr, scale):
  """Copy the elements of an array <scale> times to enlarge an image array.

  Args:
    arr: numpy array. Input array.
    scale: positive integer. Number of times to enlarge the array.

  Returns:
    numpy array.
  """
  if scale == 1:
    return arr
  tmp = np.repeat(arr, scale, axis=1)  # repeat the columns
  return np.repeat(tmp, scale, axis=0)  # repeat the rows


def _get_image_type_from_array(arr):
  """Find image type based on array dimensions.

  Raises error on invalid image dimensions.
  Args:
    arr: numpy array. Input array.

  Returns:
    str. "RGB" or "L", meant for PIL.Image.fromarray.
  """
  if len(arr.shape) == 3 and arr.shape[2] == 3:
    # 8-bit x 3 colors
    return "RGB"
  elif len(arr.shape) == 2:
    # 8-bit, gray-scale
    return "L"
  else:
    raise ValueError(
        "Input array must have either 2 dimensions or 3 dimensions where the "
        "third dimension has 3 channels. i.e. arr.shape is (x,y) or (x,y,3). "
        "Found shape {}.".format(arr.shape))


def scale_array_for_image(arr, vmin=None, vmax=None, scale=None):
  """Adjust an array to prepare it for saving to an image.

  Re-scale numbers in the input array to go from 0 to 255 to adapt them for a
  PNG image, and scale the image up to a nice size for convenience.

  Args:
    arr: numpy array. Should be 2-dimensional or 3-dimensional where the third
      dimension has 3 channels.
    vmin: number (float or int). Minimum data value, which will correspond to
      black in greyscale or lack of each color in RGB images. Default None takes
      the minimum of the data from arr.
    vmax: number (float or int). Maximum data value, which will correspond to
      white in greyscale or full presence of each color in RGB images. Default
      None takes the max of the data from arr.
    scale: integer. Number of pixels wide and tall to show each cell in the
      array. This sizes up the image while keeping exactly the same number of
      pixels for every cell in the array, preserving resolution and preventing
      any interpolation or overlapping of pixels. Default None adapts to the
      size of the image to multiply it up until a limit of 500 pixels, a
      convenient size for use in notebooks. If saving to a file for automated
      processing, scale=1 is recommended to keep output files small and simple
      while still retaining all the information content.

  Returns:
    (modified numpy array, image_mode)
  """
  image_mode = _get_image_type_from_array(arr)

  if scale is None:
    scale = max(1, int(500 / max(arr.shape)))

  if vmin is None:
    vmin = np.min(arr)
  if vmax is None:
    vmax = np.max(arr)

  # In cases where all elements are the same, fix the vmax so that even though
  # the whole image will be black, the user can at least see the shape
  if vmin == vmax:
    vmax = vmin + 1

  scaled = adjust_colors_for_png(arr, vmin=vmin, vmax=vmax)
  scaled = enlarge_image_array(scaled, scale=scale)
  return scaled, image_mode


def save_to_png(arr, path=None, image_mode=None, show=True):
  """Make a PNG and show it from a numpy array of dtype=np.uint8.

  Args:
    arr: numpy array. Input array to save.
    path: str. file path at which to save the image.
    image_mode: "RGB" or "L". Leave as default=None to choose based on image
      dimensions.
    show: bool. Whether to display the image using IPython (for notebooks).

  Returns:
    None.
  """
  if image_mode is None:
    image_mode = _get_image_type_from_array(arr)

  img = Image.fromarray(arr, mode=image_mode)

  # Saving to a temporary file is needed even when showing in a notebook
  if path is None:
    path = "/tmp/tmp.png"
  img.save(path)

  # Show image (great for notebooks)
  if show:
    display.display(display.Image(path))


def array_to_png(arr, path=None, show=True, vmin=None, vmax=None, scale=None):
  """Save an array as a PNG image with PIL and show it.

  Args:
    arr: numpy array. Should be 2-dimensional or 3-dimensional where the third
      dimension has 3 channels.
    path: str. Path for the image output. Default is /tmp/tmp.png for quickly
      showing the image in a notebook.
    show: bool. Whether to show the image using IPython utilities, only works in
      notebooks.
    vmin: number. Minimum data value, which will correspond to black in
      greyscale or lack of each color in RGB images. Default None takes the
      minimum of the data from arr.
    vmax: number. Maximum data value, which will correspond to white in
      greyscale or full presence of each color in RGB images. Default None takes
      the max of the data from arr.
    scale: integer. Number of pixels wide and tall to show each cell in the
      array. This sizes up the image while keeping exactly the same number of
      pixels for every cell in the array, preserving resolution and preventing
      any interpolation or overlapping of pixels. Default None adapts to the
      size of the image to multiply it up until a limit of 500 pixels, a
      convenient size for use in notebooks. If saving to a file for automated
      processing, scale=1 is recommended to keep output files small and simple
      while still retaining all the information content.

  Returns:
    None. Saves an image at path and optionally shows it with IPython.display.
  """
  scaled, image_mode = scale_array_for_image(
      arr, vmin=vmin, vmax=vmax, scale=scale)
  save_to_png(scaled, path=path, show=show, image_mode=image_mode)


def draw_deepvariant_pileup(example=None,
                            channels=None,
                            composite_type=None,
                            path=None,
                            show=True,
                            scale=None):
  """Quick utility for showing a pileup example as channels or RGB.

  Args:
    example: A tensorflow Example containing image/encoded and image/shape
      features. Will be parsed through channels_from_example. Ignored if
      channels are provided directly.
    channels: list of 2D arrays containing the data to draw.
    composite_type: str or None. Method for combining channels. One of
      [None,"RGB"].
    path: str. Output file path for saving as an image. If None, just show plot.
    show: bool. Whether to display the image for ipython notebooks. Set to False
      to prevent extra output when running in bulk.
    scale: integer. Multiplier to enlarge the image. Default: None, which will
      set it automatically for a human-readable size. Set to 1 for no scaling.

  Returns:
    None.
  """
  if example and not channels:
    channels = channels_from_example(example)
  elif not channels:
    raise ValueError("Either example OR channels must be specified.")

  if composite_type is None:
    img_array = np.concatenate(channels, axis=1)
  elif composite_type == "RGB":
    img_array = convert_6_channels_to_rgb(channels)
  else:
    raise ValueError(
        "Unrecognized composite_type: {}. Must be None or 'RGB'".format(
            composite_type))

  array_to_png(img_array, path=path, show=show, scale=scale)


def variant_from_example(example):
  """Extract Variant object from the 'variant/encoded' feature of an Example.

  Args:
    example: a DeepVariant-style make_examples output example.

  Returns:
    A Nucleus Variant.
  """
  features = example.features.feature
  var_string = features["variant/encoded"].bytes_list.value[0]
  return variants_pb2.Variant.FromString(var_string)


def locus_id_from_variant(variant):
  """Create a locus ID of form "chr:pos_ref" from a Variant object.

  Args:
    variant: a nucleus variant.

  Returns:
    str.
  """
  return "{}:{}_{}".format(variant.reference_name, variant.start,
                           variant.reference_bases)


def alt_allele_indices_from_example(example):
  """Extract indices of the particular alt allele(s) the example represents.

  Args:
    example: a DeepVariant make_examples output example.

  Returns:
    list of indices.
  """
  features = example.features.feature
  val = features["alt_allele_indices/encoded"].bytes_list.value[0]
  # Extract the encoded proto into unsigned integers and convert to regular ints
  mapped = [int(x) for x in np.frombuffer(val, dtype=np.uint8)]
  # Format is [<field id + type>, <number of elements in array>, ...<array>]
  # Extract the array only, leaving out the metadata
  return mapped[2:]


def alt_bases_from_indices(alt_allele_indices, alternate_bases):
  """Get alt allele bases based on their indices.

  e.g. one alt allele: [0], ["C"] => "C"
  or with two alt alleles: [0,2], ["C", "TT", "A"] => "C-A"

  Args:
    alt_allele_indices: list of integers. Indices of the alt alleles for a
      particular example.
    alternate_bases: list of strings. All alternate alleles for the variant.

  Returns:
    str. Alt allele(s) at the indices, joined by '-' if more than 1.
  """
  alleles = [alternate_bases[i] for i in alt_allele_indices]
  # avoiding '/' to support use in file paths
  return "-".join(alleles)


def alt_from_example(example):
  """Get alt allele(s) from a DeepVariant example.

  Args:
    example: a DeepVariant make_examples output example

  Returns:
    str. The bases of the alt alleles, joined by a -.
  """
  variant = variant_from_example(example)
  indices = alt_allele_indices_from_example(example)
  return alt_bases_from_indices(indices, variant.alternate_bases)


def locus_id_with_alt(example):
  """Get complete locus ID from a DeepVariant example.

  Args:
    example: a DeepVariant make_examples output example

  Returns:
    str in the form "chr:pos_ref_alt.
  """
  variant = variant_from_example(example)
  locus_id = locus_id_from_variant(variant)
  alt = alt_from_example(example)
  return "{}_{}".format(locus_id, alt)
