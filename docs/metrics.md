# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 102m
call_variants                    | 177m (206m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 78m
total                            | 357m = 6 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501484   | 3017     | 1369     | 0.99402       | 0.997387         | 0.995701        |
| SNP   | 3306946  | 20550    | 6060     | 0.993824      | 0.998172         | 0.995993        |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 8m
call_variants                    | 2m
postprocess_variants (with gVCF) | 1m
total                            | 11m

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1023     | 28       | 17       | 0.973359      | 0.983947         | 0.978624        |
| SNP   | 24964    | 315      | 167      | 0.987539      | 0.993355         | 0.990439        |


## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 120
call_variants                    | 174m (155m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 60m
total                            | 354m = 5.0 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

(The input BAM is phased already and DeepVariant was run with
`--use_hp_information=true`.)

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501999   | 2502     | 2718     | 0.995041      | 0.994828         | 0.994934        |
| SNP   | 3323838  | 3657     | 2143     | 0.998901      | 0.999356         | 0.999129        |

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 160m
call_variants                    | 188m (171m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 61m
total                            | 409 m = 6.8 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503173   | 1328     | 2147     | 0.997368      | 0.995953         | 0.996660        |
| SNP   | 3323668  | 3827     | 1921     | 0.99885       | 0.999423         | 0.999136        |

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

<a name="vfootnote1">*</a>: To use OpenVINO on Intel CPUs, run with
`--call_variants_extra_args "use_openvino=true"` with the Docker one-step
command. Also see https://github.com/google/deepvariant/pull/363 for more
details.

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
