# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | ~450m
call_variants for HG002          | ~365m
call_variants for HG003          | ~365m
call_variants for HG004          | ~365m
postprocess_variants (parallel)  | ~85m
total                            | ~1630m = 27.2 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 523008   | 2461     | 808      | 0.995317      | 0.99852          | 0.996916        |
| SNP   | 3347844  | 17283    | 3088     | 0.994864      | 0.999079         | 0.996967        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501723   | 2778     | 1218     | 0.994494      | 0.997676         | 0.996082        |
| SNP   | 3308633  | 18863    | 3376     | 0.994331      | 0.998981         | 0.996651        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 507730   | 2789     | 1130     | 0.994537      | 0.997871         | 0.996201        |
| SNP   | 3326956  | 19654    | 3229     | 0.994127      | 0.999031         | 0.996573        |

## PacBio (HiFi)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | ~840m
call_variants for HG002          | ~325m
call_variants for HG003          | ~325m
call_variants for HG004          | ~325m
postprocess_variants (parallel)  | ~80m
total                            | ~1895m = 31.6 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth). HG002/HG003/HG004 trio is used for training a PacBio model, only chr20,
chr21, and chr22 are held out. Accuracy metrics are reported for the whole
genome for consistency.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 524046   | 1423     | 1929     | 0.997292      | 0.996494         | 0.996893        |
| SNP   | 3361595  | 3532     | 1044     | 0.99895       | 0.99969          | 0.99932         |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502625   | 1876     | 2343     | 0.996281      | 0.99556          | 0.995921        |
| SNP   | 3323671  | 3824     | 1522     | 0.998851      | 0.999543         | 0.999197        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 508515   | 2004     | 2470     | 0.996075      | 0.995379         | 0.995727        |
| SNP   | 3342962  | 3648     | 1368     | 0.99891       | 0.999591         | 0.99925         |

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | --------------
make_examples                    | ~18m
call_variants for HG002          | ~6m
call_variants for HG003          | ~6m
call_variants for HG004          | ~6m
postprocess_variants (parallel)  | ~1m
total                            | ~37m

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth).

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1058     | 32       | 13       | 0.970642      | 0.988095         | 0.979291        |
| SNP   | 25165    | 249      | 28       | 0.990202      | 0.998889         | 0.994526        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1020     | 31       | 12       | 0.970504      | 0.988615         | 0.979476        |
| SNP   | 25036    | 243      | 28       | 0.990387      | 0.998883         | 0.994617        |


#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1056     | 26       | 19       | 0.97597       | 0.982852         | 0.979399        |
| SNP   | 24945    | 232      | 17       | 0.990785      | 0.999319         | 0.995034        |


## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration. For more scalable execution,
see the [External Solutions] section.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_deeptrio.sh

# WGS
bash inference_deeptrio.sh --model_preset WGS

# WES
bash inference_deeptrio.sh --model_preset WES

# PacBio
bash inference_deeptrio.sh --model_preset PACBIO

```

Runtime metrics are taken from the resulting log after each stage of DeepTrio,
and the accuracy metrics come from the hap.py summary.csv output file.

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
