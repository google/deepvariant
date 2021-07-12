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
| INDEL | 523008   | 2461     | 808      | 0.995317      | 0.998520         | 0.996916        |
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
make_examples                    | ~850m
call_variants for HG002          | ~320m
call_variants for HG003          | ~320m
call_variants for HG004          | ~320m
postprocess_variants (parallel)  | ~70m
total                            | ~1880m = 31.3 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth). HG002/HG003/HG004 trio is used for training a PacBio model, only chr20,
chr21, and chr22 are held out. Accuracy metrics are reported for the whole
genome for consistency.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 524382   | 1087     | 1518     | 0.997931      | 0.997243         | 0.997587        |
| SNP   | 3361791  | 3336     | 719      | 0.999009      | 0.999786         | 0.999397        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502854   | 1647     | 2205     | 0.996735      | 0.995824         | 0.996279        |
| SNP   | 3323866  | 3629     | 1205     | 0.998909      | 0.999638         | 0.999273        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 508752   | 1767     | 2290     | 0.996539      | 0.995717         | 0.996128        |
| SNP   | 3343130  | 3480     | 1104     | 0.998960      | 0.999670         | 0.999315        |

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
