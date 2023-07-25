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

"""Utility functions for visualization and inspection of pileup examples.

Visualization and inspection utility functions enable showing image-like array
data including those used in DeepVariant.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import enum
import math
from typing import List, NamedTuple, Tuple

from etils import epath
from IPython import display
import numpy as np
from PIL import Image
from PIL import ImageDraw

from third_party.nucleus.protos import variants_pb2


DEEPVARIANT_CHANNEL_NAMES = [
    'read base', 'base quality', 'mapping quality', 'strand',
    'read supports variant', 'base differs from ref', 'haplotype tag',
    'alternate allele 1', 'alternate allele 2'
]


class Diff(enum.Enum):
  FEW_DIFFS = 1
  MANY_DIFFS = 2
  NEARBY_VARIANTS = 3


class BaseQuality(enum.Enum):
  GOOD = 1
  BAD = 2


class MappingQuality(enum.Enum):
  GOOD = 1
  BAD = 2


class StrandBias(enum.Enum):
  GOOD = 1
  BIASED = 2


class ReadSupport(enum.Enum):
  ALL = 1
  HALF = 2
  LOW = 3


PileupCuration = NamedTuple('PileupCuration',
                            [('base_quality', BaseQuality),
                             ('mapping_quality', MappingQuality),
                             ('strand_bias', StrandBias),
                             ('diff_category', Diff),
                             ('read_support', ReadSupport)])


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
  img = features['image/encoded'].bytes_list.value[0]
  shape = features['image/shape'].int64_list.value[0:3]
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


def scale_colors_for_png(arr, vmin=0, vmax=255):
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
    raise ValueError('vmin must be non-zero and higher than vmin.')

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
    return 'RGB'
  elif len(arr.shape) == 2:
    # 8-bit, gray-scale
    return 'L'
  else:
    raise ValueError(
        'Input array must have either 2 dimensions or 3 dimensions where the '
        'third dimension has 3 channels. i.e. arr.shape is (x,y) or (x,y,3). '
        'Found shape {}.'.format(arr.shape))


def autoscale_colors_for_png(arr, vmin=None, vmax=None):
  """Adjust an array to prepare it for saving to an image.

  Re-scale numbers in the input array to go from 0 to 255 to adapt them for a
  PNG image.

  Args:
    arr: numpy array. Should be 2-dimensional or 3-dimensional where the third
      dimension has 3 channels.
    vmin: number (float or int). Minimum data value, which will correspond to
      black in greyscale or lack of each color in RGB images. Default None takes
      the minimum of the data from arr.
    vmax: number (float or int). Maximum data value, which will correspond to
      white in greyscale or full presence of each color in RGB images. Default
      None takes the max of the data from arr.

  Returns:
    (modified numpy array, image_mode)
  """
  image_mode = _get_image_type_from_array(arr)

  if vmin is None:
    vmin = np.min(arr)
  if vmax is None:
    vmax = np.max(arr)

  # In cases where all elements are the same, fix the vmax so that even though
  # the whole image will be black, the user can at least see the shape
  if vmin == vmax:
    vmax = vmin + 1

  scaled = scale_colors_for_png(arr, vmin=vmin, vmax=vmax)
  return scaled, image_mode


def add_header(img, labels, mark_midpoints=True, header_height=20):
  """Adds labels to the image, evenly distributed across the top.

  This is primarily useful for showing the names of channels.

  Args:
    img: A PIL Image.
    labels: list of strs. Labels for segments to write across the top.
    mark_midpoints: bool. Whether to add a small vertical line marking the
      center of each segment of the image.
    header_height: int. Height of the header in pixels.

  Returns:
    A new PIL Image, taller than the original img and annotated.
  """

  # Create a taller image to make space for a header at the top.
  new_height = header_height + img.size[1]
  new_width = img.size[0]

  if img.mode == 'RGB':
    placeholder_size = (new_height, new_width, 3)
  else:
    placeholder_size = (new_height, new_width)
  placeholder = np.ones(placeholder_size, dtype=np.uint8) * 255

  # Divide the image width into segments.
  segment_width = img.size[0] / len(labels)

  # Calculate midpoints for all segments.
  midpoints = [int(segment_width * (i + 0.5)) for i in range(len(labels))]

  if mark_midpoints:
    # For each label, add a small line to mark the middle.
    for x_position in midpoints:
      placeholder[header_height - 5:header_height, x_position] = 0
      # If image has an even width, it will need 2 pixels marked as the middle.
      if segment_width % 2 == 0:
        placeholder[header_height - 5:header_height, x_position + 1] = 0

  bigger_img = Image.fromarray(placeholder, mode=img.mode)
  # Place the original image inside the taller placeholder image.
  bigger_img.paste(img, (0, header_height))

  # Add a label for each segment.
  draw = ImageDraw.Draw(bigger_img)
  for i in range(len(labels)):
    text = labels[i]
    text_width = draw.textbbox((0, 0), text, anchor='lt')[2]
    # xy refers to the left top corner of the text, so to center the text on
    # the midpoint, subtract half the text width from the midpoint position.
    x_position = int(midpoints[i] - text_width / 2)
    draw.text(xy=(x_position, 0), text=text, fill='black')
  return bigger_img


def save_to_png(arr,
                path=None,
                image_mode=None,
                show=True,
                labels=None,
                scale=None):
  """Make a PNG and show it from a numpy array of dtype=np.uint8.

  Args:
    arr: numpy array. Input array to save.
    path: str. File path at which to save the image. A .png prefix is added if
      the path does not already have one. Leave empty to save at /tmp/tmp.png,
      which is useful when only temporarily showing the image in a Colab
      notebook.
    image_mode: "RGB" or "L". Leave as default=None to choose based on image
      dimensions.
    show: bool. Whether to display the image using IPython (for notebooks).
    labels: list of str. Labels to show across the top of the image.
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
  if image_mode is None:
    image_mode = _get_image_type_from_array(arr)

  img = Image.fromarray(arr, mode=image_mode)

  if labels is not None:
    img = add_header(img, labels)

  if scale is None:
    scale = max(1, int(500 / max(arr.shape)))

  if scale != 1:
    img = img.resize((img.size[0] * scale, img.size[1] * scale))

  # Saving to a temporary file is needed even when showing in a notebook
  if path is None:
    path = '/tmp/tmp.png'
  elif not path.endswith('.png'):
    # Only PNG is supported because JPEG files are unnecessarily 3 times larger.
    path = '{}.png'.format(path)
  img.save(epath.Path(path).open('wb'), format=path.split('.')[-1])

  # Show image (great for notebooks)
  if show:
    display.display(display.Image(path))


def array_to_png(arr,
                 path=None,
                 show=True,
                 vmin=None,
                 vmax=None,
                 scale=None,
                 labels=None):
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
    labels: list of str. Labels to show across the top of the image.

  Returns:
    None. Saves an image at path and optionally shows it with IPython.display.
  """
  scaled, image_mode = autoscale_colors_for_png(arr, vmin=vmin, vmax=vmax)
  save_to_png(
      scaled,
      path=path,
      show=show,
      image_mode=image_mode,
      labels=labels,
      scale=scale)


def _deepvariant_channel_names(num_channels):
  """Get DeepVariant channel names for the given number of channels."""
  # Add additional empty labels if there are more channels than expected.
  filler_labels = [
      'channel {}'.format(i + 1)
      for i in range(len(DEEPVARIANT_CHANNEL_NAMES), num_channels)
  ]
  labels = DEEPVARIANT_CHANNEL_NAMES + filler_labels
  # Trim off any extra labels.
  return labels[0:num_channels]


def draw_deepvariant_pileup(example=None,
                            channels=None,
                            composite_type=None,
                            annotated=True,
                            labels=None,
                            path=None,
                            show=True,
                            scale=None):
  """Quick utility for showing a pileup example as channels or RGB.

  Args:
    example: A tensorflow Example containing image/encoded and image/shape
      features. Will be parsed through channels_from_example. Ignored if
      channels are provided directly. Either example OR channels is required.
    channels: list of 2D arrays containing the data to draw. Either example OR
      channels is required.
    composite_type: str or None. Method for combining channels. One of
      [None,"RGB"].
    annotated: bool. Whether to add channel labels and mark midpoints.
    labels: list of str. Which labels to add to the image. If annotated=True,
      use default channels labels for DeepVariant.
    path: str. Output file path for saving as an image. If None, just show plot.
    show: bool. Whether to display the image for ipython notebooks. Set to False
      to prevent extra output when running in bulk.
    scale: integer. Multiplier to enlarge the image. Default: None, which will
      set it automatically for a human-readable size. Set to 1 for no scaling.

  Returns:
    None. Saves an image at path and optionally shows it with IPython.display.
  """
  if example and not channels:
    channels = channels_from_example(example)
  elif not channels:
    raise ValueError('Either example OR channels must be specified.')

  if composite_type is None:
    img_array = np.concatenate(channels, axis=1)
    if annotated and labels is None:
      labels = _deepvariant_channel_names(len(channels))
  elif composite_type == 'RGB':
    img_array = convert_6_channels_to_rgb(channels)
    if annotated and labels is None:
      labels = ['']  # Creates one midpoint with no label.
  else:
    raise ValueError(
        "Unrecognized composite_type: {}. Must be None or 'RGB'".format(
            composite_type))

  array_to_png(
      img_array,
      path=path,
      show=show,
      scale=scale,
      labels=labels,
      vmin=0,
      vmax=254)


def variant_from_example(example):
  """Extract Variant object from the 'variant/encoded' feature of an Example.

  Args:
    example: a DeepVariant-style make_examples output example.

  Returns:
    A Nucleus Variant.
  """
  features = example.features.feature
  var_string = features['variant/encoded'].bytes_list.value[0]
  return variants_pb2.Variant.FromString(var_string)


def locus_id_from_variant(variant):
  """Create a locus ID of form "chr:pos_ref" from a Variant object.

  Args:
    variant: a nucleus variant.

  Returns:
    str.
  """
  return '{}:{}_{}'.format(variant.reference_name, variant.start,
                           variant.reference_bases)


def alt_allele_indices_from_example(example):
  """Extract indices of the particular alt allele(s) the example represents.

  Args:
    example: a DeepVariant make_examples output example.

  Returns:
    list of indices.
  """
  features = example.features.feature
  val = features['alt_allele_indices/encoded'].bytes_list.value[0]
  # Extract the encoded proto into unsigned integers and convert to regular ints
  mapped = [int(x) for x in np.frombuffer(val, dtype=np.uint8)]
  # Format is [<field id + type>, <number of elements in array>, ...<array>].
  # Extract the array only, leaving out the metadata.
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
  # Avoiding '/' to support use in file paths.
  return '-'.join(alleles)


def alt_from_example(example):
  """Get alt allele(s) from a DeepVariant example.

  Args:
    example: a DeepVariant make_examples output example.

  Returns:
    str. The bases of the alt alleles, joined by a -.
  """
  variant = variant_from_example(example)
  indices = alt_allele_indices_from_example(example)
  return alt_bases_from_indices(indices, variant.alternate_bases)


def locus_id_with_alt(example):
  """Get complete locus ID from a DeepVariant example.

  Args:
    example: a DeepVariant make_examples output example.

  Returns:
    str in the form "chr:pos_ref_alt.
  """
  variant = variant_from_example(example)
  locus_id = locus_id_from_variant(variant)
  alt = alt_from_example(example)
  return '{}_{}'.format(locus_id, alt)


def label_from_example(example):
  """Get the "label" from an example.

  Args:
    example: a DeepVariant make_examples output example.

  Returns:
    integer (0, 1, or 2 for regular DeepVariant examples) or None if the
        example has no label.
  """
  val = example.features.feature['label'].int64_list.value
  if val:
    return int(val[0])
  else:
    return None


def remove_ref_band(arr: np.ndarray,
                    num_top_rows_to_skip: int = 5) -> np.ndarray:
  """Removes the reference rows at the top of a pileup image array."""
  assert len(arr.shape) == 2
  assert arr.shape[0] > num_top_rows_to_skip
  return arr[num_top_rows_to_skip:, :]


def fraction_low_base_quality(channels: List[np.ndarray],
                              threshold: int = 127) -> float:
  """Gets fraction of bases that have low base quality scores in a pileup.

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[1], the base quality channel.
      threshold: Bases qualities below this will be considered low quality. The
        default is 127 because this is half of the max (254).

  Returns:
      The fraction of bases with base quality below the threshold.
  """
  basequal_channel = remove_ref_band(channels[1])
  non_zero_values = basequal_channel[basequal_channel > 0]

  num_non_zero = non_zero_values.shape[0]
  if num_non_zero == 0:
    return 0.0
  return sum((non_zero_values < threshold) * 1.0) / num_non_zero


def fraction_reads_with_low_mapq(channels: List[np.ndarray],
                                 threshold: int = 127) -> float:
  """Gets fraction of reads that have low mapping quality scores in pileup.

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[2], the mapping quality channel.
      threshold: int. Default is 127 because this is half of the max (254).

  Returns:
      The fraction of bases with mapping quality below the threshold.
  """
  mapq_channel = remove_ref_band(channels[2])
  # Get max value of each row, aka each read.
  max_row_values = np.amax(mapq_channel, axis=1)

  non_zero_values = max_row_values[max_row_values > 0]
  num_non_zero = non_zero_values.shape[0]
  if num_non_zero == 0:
    return 0.0
  return sum((non_zero_values < threshold) * 1.0) / num_non_zero


def fraction_read_support(channels: List[np.ndarray]) -> float:
  """Gets fraction of reads that support the variant.

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[4], the 'read supports variant' channel.

  Returns:
      Fraction of reads supporting the alternate allele(s), ranging from [0, 1].
  """
  support_channel = remove_ref_band(channels[4])
  max_row_values = np.amax(support_channel, axis=1)

  non_zero_values = max_row_values[max_row_values > 0]
  num_non_zero = non_zero_values.shape[0]
  if num_non_zero == 0:
    return 0.0
  return sum(non_zero_values == 254) * 1.0 / num_non_zero


def describe_read_support(channels: List[np.ndarray]) -> ReadSupport:
  """Calculates read support and describes it categorically.

  Computes read support as a fraction and returns a convenient descriptive term
  according to the following thresholds: LOW is [0, 0.3], HALF is (0.3, 0.8],
  and ALL is (0.8, 1].

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[4], the 'read supports variant' channel.

  Returns:
      A ReadSupport value.
  """
  fraction_support = fraction_read_support(channels)
  if fraction_support > 0.8:
    return ReadSupport.ALL
  elif fraction_support > 0.3:
    return ReadSupport.HALF
  else:
    return ReadSupport.LOW


def binomial_test(k: int, n: int) -> float:
  """Calculates a two-tailed binomial test with p=0.5, without scipy.

  Since the expected probability is 0.5, it simplifies a few things:
  1) (0.5**x)*(0.5**(n-x)) = (0.5**n)
  2) A two-tailed test is simply doubling when p = 0.5.
  Scipy is much larger than Nucleus, so this avoids adding it as a dependency.

  Args:
    k: Number of "successes", in this case, the number of supporting reads.
    n: Number of "trials", in this case, the total number of reads.

  Returns:
    The p-value for the binomial test.
  """
  if not k <= n:
    raise ValueError('k must be <= n')
  if k == n / 2:
    return 1.0
  sum_of_ps = 0

  # With p=0.5, the distribution is symmetric, allowing this simplification:
  k = min(k, n - k)
  # Add up all the exact probabilities for each scenario more extreme than k.
  for x in range(0, k + 1):
    # After python 3.8, the following line can be replaced using math.comb.
    n_choose_x = math.factorial(n) / math.factorial(x) / math.factorial(n - x)
    p_for_i = n_choose_x * (0.5**n)
    sum_of_ps += p_for_i
  return sum_of_ps * 2  # Doubling because it's a two-tailed test.


def pvalue_for_strand_bias(channels: List[np.ndarray]) -> float:
  """Calculates a rough p-value for strand bias in pileup.

  Using the strand and read-supports-variant channels, compares the numbers of
  forward and reverse reads among the supporting reads and returns a p-value
  using a two-tailed binomial test.

  Args:
      channels: List of DeepVariant channels. Uses channels[3] (strand) and
        channels[4] (read support).

  Returns:
      P-value for whether the supporting reads show strand bias.
  """
  strand = remove_ref_band(channels[3])
  forward_strand = strand == 240
  reverse_strand = strand == 70

  read_support = remove_ref_band(channels[4])
  read_support = (read_support == 254) * 1.0
  forward_support = read_support * forward_strand
  reverse_support = read_support * reverse_strand

  forward_supporting = int(sum(np.amax(forward_support, axis=1)))
  reverse_supporting = int(sum(np.amax(reverse_support, axis=1)))

  return binomial_test(
      k=forward_supporting, n=forward_supporting + reverse_supporting)


def analyze_diff_and_nearby_variants(
    channels: List[np.ndarray]) -> Tuple[float, int]:
  """Analyzes which differences belong to nearby variants and which do not.

  This attempts to identify putative nearby variants from the pileup image
  alone, and then excludes these columns of the pileup to calculate the
  remaining fraction of differences that may be attributed to sequencing errors.

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[5], the 'differs from ref' channel.

  Returns:
      Two outputs: diff fraction, number of likely nearby variants.
  """
  diff_channel = remove_ref_band(channels[5])

  # Count the number of diff pixels per column.
  column_diffs = np.sum(diff_channel == 254, axis=0)
  # Count number of differences per base position.
  column_read_count = np.sum(diff_channel != 0, axis=0)
  # Divide to get the fraction of reads showing a diff at each base (column).
  # Adding 1 here avoids dividing by zero (exact fraction here is not vital).
  fraction = column_diffs * 1.0 / (column_read_count + 1)

  # Columns with more differences could be variants.
  nearby_variant_columns = (fraction > 0.1) * (column_diffs > 4) * 1
  num_potential_nearby_variants = sum(nearby_variant_columns)

  # Exclude potential variants when calculating error fraction.
  nearby_variant_mask = np.array([nearby_variant_columns] *
                                 diff_channel.shape[0])
  mask_to_remove_nearby_variants = 1 - nearby_variant_mask
  non_variant_diffs = (diff_channel == 254) * mask_to_remove_nearby_variants

  # Calculate differences as fraction of the total number of read bases.
  total_read_area = np.sum((diff_channel != 0))
  diff_fraction = 0 if total_read_area == 0 else np.sum(
      non_variant_diffs) / total_read_area
  return diff_fraction, num_potential_nearby_variants


def describe_diff(channels: List[np.ndarray],
                  diff_fraction_threshold: float = 0.01) -> Diff:
  """Describes a pileup image by its diff channel, including nearby variants.

  Returns Diff.MANY_DIFFS if the fraction of differences outside potential
  nearby variants is above the diff_fraction_threshold, which is usually
  indicative of sequencing errors. Otherwise return Diff.NEARBY_VARIANTS if
  there are five or more of these, or Diff.FEW_DIFFS if neither of these
  special cases apply.

  Args:
      channels: A list of channels of a DeepVariant pileup image. This only uses
        channels[5], the 'differs from ref' channel.
      diff_fraction_threshold: Fraction of total bases of all reads that can
        differ, above which the pileup will be designated as 'many_diffs'.
        Differences that appear due to nearby variants (neater columns) do not
        count towards this threshold. The default is set by visual curation of
        Illumina reads, so it may be necessary to increase this for higher-error
        sequencing types.

  Returns:
      One Diff value.
  """
  diff_fraction, nearby_variants = analyze_diff_and_nearby_variants(channels)
  # Thresholds were chosen by visual experimentation, i.e. human curation.
  if diff_fraction > diff_fraction_threshold:
    return Diff.MANY_DIFFS
  elif nearby_variants >= 5:
    return Diff.NEARBY_VARIANTS
  else:
    return Diff.FEW_DIFFS


def curate_pileup(channels: List[np.ndarray]) -> PileupCuration:
  """Runs all automated curation functions and outputs categorical tags.

  The following values are possible for each descriptor:
  * base_quality: GOOD (>5% low quality) or BAD.
  * mapping_quality: GOOD (<5% low quality) or BAD.
  * strand_biased: BIASED (p-value < 0.05) or GOOD.
  * diff_category: MANY_DIFFS (>1% differences), NEARBY_VARIANTS (5+ variants),
  or FEW_DIFFS otherwise.
  * read_support: LOW (<=30%), HALF (30-80%), ALL (>80%).

  The thresholds were all set by trying to match human curation.

  Args:
      channels: A list of DeepVariant channels.

  Returns:
      A PileupCuration NamedTuple.
  """

  return PileupCuration(
      base_quality=BaseQuality.GOOD
      if fraction_low_base_quality(channels) < 0.05 else BaseQuality.BAD,
      mapping_quality=MappingQuality.GOOD
      if fraction_reads_with_low_mapq(channels) < 0.05 else MappingQuality.BAD,
      strand_bias=StrandBias.BIASED
      if pvalue_for_strand_bias(channels) < 0.05 else StrandBias.GOOD,
      diff_category=describe_diff(channels),
      read_support=describe_read_support(channels))
