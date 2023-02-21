# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~98m
call_variants                    | ~174m
postprocess_variants (with gVCF) | ~46m
total                            | ~318m = ~5.3 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501712   | 2789     | 1175     | 0.994472      | 0.997758         | 0.996112        |
| SNP   | 3306963  | 20533    | 3816     | 0.993829      | 0.998848         | 0.996332        |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~6m
call_variants                    | ~1m
postprocess_variants (with gVCF) | ~1m
total                            | ~8m

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1019     | 32       | 10       | 0.969553      | 0.990467         | 0.979898        |
| SNP   | 24981    | 298      | 49       | 0.988212      | 0.998043         | 0.993103        |

## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~153m
call_variants                    | ~196m
postprocess_variants (with gVCF) | ~53m
total                            | ~402m = ~6.7 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501620   | 2881     | 2803     | 0.994289      | 0.994664         | 0.994477        |
| SNP   | 3324720  | 2775     | 1887     | 0.999166      | 0.999433         | 0.999300        |

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~144m
call_variants                    | ~173m
postprocess_variants (with gVCF) | ~39m
total                            | ~356m = ~5.93 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503347   | 1154     | 2003     | 0.997713      | 0.996225         | 0.996968        |
| SNP   | 3323945  | 3550     | 1535     | 0.998933      | 0.999539         | 0.999236        |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.4/scripts/inference_deepvariant.sh

# WGS
bash inference_deepvariant.sh --model_preset WGS

# WES
bash inference_deepvariant.sh --model_preset WES

# PacBio
bash inference_deepvariant.sh --model_preset PACBIO

# Hybrid
bash inference_deepvariant.sh --model_preset HYBRID_PACBIO_ILLUMINA
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant. The runtime numbers reported above are the average of 5 runs each.
The accuracy metrics come from the hap.py summary.csv output file.
The runs are deterministic so all 5 runs produced the same output.

[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
