# Copyright 2024 Google LLC.
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
"""Module containing helper functions for running small model training."""

import os
import time
from typing import Any, Sequence, Tuple, Union

from absl import logging
import matplotlib.pyplot as plt
import ml_collections
import numpy as np
import pandas as pd
import pysam
from sklearn import ensemble
import tensorflow as tf

from deepvariant.small_model import make_small_model_examples
from io import fileio
from third_party.nucleus.util import genomics_math


def _is_snp(variant: pysam.VariantRecord) -> bool:
  """Determines if the given VariantRecord is a SNP."""
  return len(variant.ref) == 1 and all(len(x) == 1 for x in variant.alts)


def _get_genotype_class(genotype: Tuple[int, int]) -> int:
  """Safely tries to encode the genotype class for the given genotype."""
  try:
    return make_small_model_examples.ENCODING_BY_GENOTYPE[genotype]
  except KeyError:
    return -1


def _format_variant_to_dict(variant: pysam.VariantRecord) -> dict[str, Any]:
  """Formats the given VariantRecord to a dict."""
  variant_dict = {
      "contig": variant.contig,
      "start": variant.start,
      "stop": variant.stop,
      "qual": variant.qual,
      "ref": variant.ref,
      "alts": variant.alts,
      "id": f"{variant.contig}_{variant.start}_{variant.stop}",
  }
  for index, sample_name in enumerate(variant.samples):
    sample = variant.samples[index]
    for key, value in sample.items():
      variant_dict[f"{sample_name}_{key}"] = value
    if None in sample["GT"]:
      genotype = None
    else:
      gt_0, gt_1 = sample["GT"]
      genotype = (
          variant.alleles[gt_0],
          variant.alleles[gt_1],
      )
    variant_dict[f"{sample_name}_genotype"] = genotype
    variant_dict[f"{sample_name}_gt"] = _get_genotype_class(sample["GT"])
  return variant_dict


def load_vcf_to_df(
    variant_file_path: str, regions: list[str] | None = None
) -> pd.DataFrame:
  """Loads a VCF to dataframe, caching the result to speed up repeated runs."""
  cache_vcf_csv = variant_file_path + "_cache.csv"
  if fileio.Exists(cache_vcf_csv):
    logging.info("Loading VCF from cache")
    return pd.read_csv(cache_vcf_csv, index_col=0)

  logging.info("Loading VCF to DataFrame")
  variants = []
  for variant in pysam.VariantFile(variant_file_path):
    if regions and variant.contig not in regions:
      continue
    if not _is_snp(variant):
      continue
    variants.append(_format_variant_to_dict(variant))
  df = pd.DataFrame(variants)
  # Index by the `id` value
  df.set_index("id", inplace=True)
  df.to_csv(cache_vcf_csv)
  return df


def load_tsv(path: str) -> pd.DataFrame:
  """Loads a small model example TSV file."""
  try:
    return pd.read_csv(path, sep="\t", header=0)
  except pd.errors.EmptyDataError:
    return pd.DataFrame()


def load_files_from_directory_tree(
    tsv_directory: str, max_num_records: Union[None, int] = None
) -> pd.DataFrame:
  """Recursively loads all *small_model.tsv files in the directory."""
  dfs = []
  current_num_records = 0
  for directory, _, filenames in fileio.Walk(tsv_directory):
    for filename in filenames:
      if filename.endswith("small_model.tsv"):
        logging.info("Loading %s", os.path.join(directory, filename))
        with fileio.Open(os.path.join(directory, filename), "r") as f:
          df = load_tsv(f)
          current_num_records += len(df)
          dfs.append(df)
          if max_num_records and current_num_records >= max_num_records:
            break
  return pd.concat(dfs).reset_index(drop=True)


def filter_out_mislabels(vcf: pd.DataFrame, df: pd.DataFrame) -> pd.DataFrame:
  """Given a truth set, filter out any candidate example that is mislabeled."""
  logging.info("Calculating mislabels in dataframe.")
  gt_column = [c for c in vcf.columns if c.endswith("_gt")][0]

  def get_true_gt(row):
    key = f"{row.contig}_{row.start}_{row.end}"
    try:
      return vcf.loc[key][gt_column]
    except KeyError:
      # if not present, must be a no-call
      return make_small_model_examples.GenotypeEncoding.REF.value

  df["true_gt"] = df.apply(get_true_gt, axis=1)
  mislabels = df[df["true_gt"] != df["genotype"]]
  df.drop("true_gt", axis=1, inplace=True)
  return df.drop(mislabels.index)


def log_dataframe(df: pd.DataFrame, title: str = "") -> None:
  """Inserts a crucial newline to keep df alignment for logging."""
  logging.info("\n".join(map(str, (title, df))))


def sanitize_mislabels(
    truth_vcf: str, candidates: pd.DataFrame
) -> pd.DataFrame:
  """Sanitizes mislabels in the candidates dataframe."""
  logging.info("Attempting to filter mislabeled examples.")
  vcf = load_vcf_to_df(truth_vcf)
  logging.info("Loaded VCF")
  log_dataframe(vcf, "Truth VCF")
  sanitized_candidates = filter_out_mislabels(vcf, candidates)
  num_filtered_examples = len(candidates) - len(sanitized_candidates)
  logging.info("Filtered out %s examples", num_filtered_examples)
  return sanitized_candidates


def load_training_data(
    tsv_directory: str, max_num_records: Union[None, int] = None
) -> pd.DataFrame:
  """Loads data from a set of TSV files and truncates to size."""
  logging.info("Loading candidates from TSV")
  candidates = load_files_from_directory_tree(tsv_directory, max_num_records)
  logging.info("Num candidates: %s", len(candidates))
  df = candidates
  if max_num_records and len(candidates) > max_num_records:
    df = candidates.sample(n=max_num_records)
  return df


def split_training_data(
    df: pd.DataFrame, test_fraction: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
  """Splits the dataframe into training and testing data."""
  df_train = df.sample(frac=(1 - test_fraction))
  df_test = df[~df.index.isin(df_train.index)]
  log_dataframe(df_train, "Training data")
  log_dataframe(df_test, "Testing data")
  return df_train, df_test


def keras_mlp_model(model_params: ml_collections.ConfigDict) -> tf.keras.Model:
  """Creates a Keras MLP model."""
  model = tf.keras.Sequential()
  input_shape = len(make_small_model_examples.MODEL_FEATURES)
  hidden_layers = model_params.hidden_layer_sizes
  model.add(
      tf.keras.layers.Dense(
          hidden_layers[0],
          activation=model_params.activation,
          input_shape=(input_shape,),
      )
  )
  if len(hidden_layers) > 1:
    for layer_size in hidden_layers[1:]:
      model.add(
          tf.keras.layers.Dense(layer_size, activation=model_params.activation)
      )
  output_shape = len(make_small_model_examples.ENCODING_BY_GENOTYPE)
  model.add(tf.keras.layers.Dense(output_shape, activation="softmax"))

  model.summary()
  model.compile(
      optimizer=model_params.optimizer,
      loss="categorical_crossentropy",
      metrics=["accuracy"],
  )
  return model


def one_hot_encode_truth(truth: pd.Series) -> np.ndarray:
  """Converts a series of encoded genotype values to one-hot encoding."""
  one_hot_encoded_truth = np.zeros(shape=(len(truth), 3))
  for index, encoded_genotype in enumerate(truth):
    one_hot_encoded_truth[index][encoded_genotype - 1] = 1
  return one_hot_encoded_truth


class ModelRunner:
  """Run model training and evaluation."""

  def __init__(
      self,
      model: tf.keras.Model,
      training_df: pd.DataFrame,
      test_df: pd.DataFrame,
      epochs: int,
      batch_size: int,
      name="",
      model_features=None,
  ):
    self.model = model
    self.training_df = training_df
    self.test_df = test_df
    self.test_df_mut = test_df.copy()
    self.name = name
    self.training_time = 0
    self.model_features = model_features
    self.inference_time = 0
    self.epochs = epochs
    self.batch_size = batch_size

  def train(self) -> None:
    """Trains the model."""
    x_train = self.training_df[self.get_model_features()]
    y_train = self.training_df[make_small_model_examples.TRUTH_FEATURE.value]
    y_train_one_hot = one_hot_encode_truth(y_train)
    time_before_train = time.time()
    self.model.fit(
        x_train.values,
        y_train_one_hot,
        epochs=self.epochs,
        batch_size=self.batch_size,
    )
    self.training_time = time.time() - time_before_train

  def test(self) -> np.ndarray:
    """Runs the inference and adds predictions to test_df_mut."""
    before_predictions = time.time()
    testing_df = self.test_df[self.get_model_features()]
    probabilities = self.model.predict(testing_df.values)
    predictions = np.array([np.argmax(p) for p in probabilities])
    self.inference_time = time.time() - before_predictions
    self.test_df_mut["gt_predicted"] = predictions
    self.test_df_mut["gq"] = [
        genomics_math.ptrue_to_bounded_phred(max(p)) for p in probabilities
    ]
    return predictions

  def assess(self, predictions: np.ndarray) -> None:
    """Runs quick analysis on model performance."""
    y_test = self.test_df[make_small_model_examples.TRUTH_FEATURE.value]
    cm = tf.math.confusion_matrix(labels=y_test, predictions=predictions)
    model_metrics = [
        f"Name: {self.get_experiment_name()}",
        f"Accuracy, {self.model.metrics[1].result()}",
        f"Num Training Candidates, {len(self.training_df)}",
        f"Num Inference Candidates, {len(self.test_df)}",
        f"Training Runtime, {self.training_time:.4f}s",
        f"Inference Runtime, {self.inference_time:.4f}s",
        "Confusion Matrix:",
        str(cm),
    ]
    logging.info("\n".join(model_metrics))

  def get_experiment_name(self) -> str:
    """Returns the experiment name, based on the classifier used and an identifier."""
    return f"{self.model.__class__.__name__.lower()}{self.name}"

  def get_model_features(self) -> Sequence[str]:
    """Returns the feature list passed to the model for this runner.

    They can be explicitly passed to this class, otherwise it will default to
    all features declared in the small_model_make_examples module that are not
    truth or serve as ids.
    """
    if self.model_features is not None:
      return self.model_features
    return [f.value for f in make_small_model_examples.MODEL_FEATURES]

  def run(self) -> None:
    """Runs model training and evaluation."""
    logging.info("Running model training and eval.")
    self.train()
    predictions = self.test()
    self.assess(predictions)

  def get_mistakes(self) -> pd.DataFrame:
    """Returns all misclassifications of the model on the test set."""
    dft = self.test_df_mut
    return dft[dft["gt_predicted"] != dft["genotype"]]

  def show_mistakes_above_threshold(self, threshold: float) -> None:
    """Logs a dataframe of all mistakes above the given GQ threshold."""
    mistakes = self.get_mistakes()
    high_gq_mistakes = mistakes[mistakes["gq"] > threshold]
    log_dataframe(high_gq_mistakes, "High Confidence Mistakes")

  def calculate_feature_importance(self) -> None:
    """Calculates feature importance using RFC and MDI."""
    x_train = self.training_df[self.get_model_features()]
    y_train = self.training_df[make_small_model_examples.TRUTH_FEATURE.value]
    feature_names = x_train.columns
    forest = ensemble.RandomForestClassifier(random_state=0)
    forest.fit(x_train, y_train)
    importances = forest.feature_importances_
    std = np.std(
        [tree.feature_importances_ for tree in forest.estimators_], axis=0
    )
    forest_importances = pd.Series(importances, index=feature_names)
    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    timestamp = str(time.time()).split(".")[0]
    png_path = f"/tmp/feature_importances_{timestamp}.png"
    logging.info("Saving feature importances to %s", png_path)
    fig.savefig(png_path)
