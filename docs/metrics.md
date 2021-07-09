# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 94m
call_variants                    | 166m
postprocess_variants (with gVCF) | 70m
total                            | 330m = 5.5 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501523   | 2978     | 1207     | 0.994097      | 0.997696         | 0.995893        |
| SNP   | 3306397  | 21099    | 4556     | 0.993659      | 0.998625         | 0.996136        |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 6m
call_variants                    | 2m
postprocess_variants (with gVCF) | 1m
total                            | 9m

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1020     | 31       | 14       | 0.970504      | 0.986717         | 0.978544        |
| SNP   | 24938    | 341      | 58       | 0.986511      | 0.997680         | 0.992064        |


## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 110m
call_variants                    | 153m
postprocess_variants (with gVCF) | 62m
total                            | 325m = 5.42 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

(The input BAM is phased already and DeepVariant was run with
`--use_hp_information=true`.)

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501805   | 2696     | 2661     | 0.994656      | 0.994935         | 0.994795        |
| SNP   | 3323555  | 3940     | 1642     | 0.998816      | 0.999507         | 0.999161        |

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 154m
call_variants                    | 180m
postprocess_variants (with gVCF) | 56m
total                            | 390m = 6.5 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503228   | 1273     | 1990     | 0.997477      | 0.996249         | 0.996863        |
| SNP   | 3323696  | 3799     | 1710     | 0.998858      | 0.999486         | 0.999172        |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration. For more scalable execution
of DeepVariant see the [External Solutions] section.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

redacted

```
# WGS
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_wgs.sh
bash inference_wgs.sh

# WES
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_wes.sh
bash inference_wes.sh

# PacBio
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_pacbio.sh
bash inference_pacbio.sh

# Hybrid
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_hybrid_pacbio_illumina.sh
bash inference_hybrid_pacbio_illumina.sh
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant, and the accuracy metrics come from the hap.py summary.csv output
file.

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
