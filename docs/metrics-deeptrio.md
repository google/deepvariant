# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | ~410m
call_variants for HG002          | ~350m
call_variants for HG003          | ~350m
call_variants for HG004          | ~350m
postprocess_variants (parallel)  | ~90m
total                            | ~1550m = ~25.8 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 523025   | 2444     | 799      | 0.995349      | 0.998537         | 0.99694         |
| SNP   | 3347837  | 17290    | 3071     | 0.994862      | 0.999084         | 0.996969        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501742   | 2759     | 1210     | 0.994531      | 0.997691         | 0.996109        |
| SNP   | 3308645  | 18851    | 3367     | 0.994335      | 0.998984         | 0.996654        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 507704   | 2815     | 1148     | 0.994486      | 0.997837         | 0.996159        |
| SNP   | 3326963  | 19647    | 3169     | 0.994129      | 0.999049         | 0.996583        |

## PacBio (HiFi)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | ~690m
call_variants for HG002          | ~300m
call_variants for HG003          | ~300m
call_variants for HG004          | ~300m
postprocess_variants (parallel)  | ~70m
total                            | ~1660m = ~27.67 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth). HG002/HG003/HG004 trio is used for training a PacBio model, only chr20,
chr21, and chr22 are held out. Accuracy metrics are reported for the whole
genome for consistency.

#### HG002:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 524377   | 1092     | 1510     | 0.997922      | 0.997257         | 0.997589        |
| SNP   | 3361785  | 3342     | 727      | 0.999007      | 0.999784         | 0.999395        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502844   | 1657     | 2209     | 0.996716      | 0.995816         | 0.996266        |
| SNP   | 3323870  | 3625     | 1198     | 0.998911      | 0.99964          | 0.999275        |

#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 508733   | 1786     | 2306     | 0.996502      | 0.995687         | 0.996094        |
| SNP   | 3343146  | 3464     | 1101     | 0.998965      | 0.999671         | 0.999318        |

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
| INDEL | 1059     | 31       | 11       | 0.97156       | 0.989918         | 0.980653        |
| SNP   | 25163    | 251      | 31       | 0.990124      | 0.99877          | 0.994428        |

#### HG003:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1018     | 33       | 14       | 0.968601      | 0.986679         | 0.977557        |
| SNP   | 25034    | 245      | 27       | 0.990308      | 0.998923         | 0.994597        |


#### HG004:

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1055     | 27       | 17       | 0.975046      | 0.984601         | 0.979801        |
| SNP   | 24947    | 230      | 17       | 0.990865      | 0.999319         | 0.995074        |


## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration. For more scalable execution,
see the [External Solutions] section.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.2/scripts/inference_deeptrio.sh

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
