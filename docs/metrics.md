# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 97m
call_variants                    | 249m (211m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 88m
total                            | 434m = 7.2 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel |  501470 |  3031 | 1380 | 0.993992 | 0.997367  | 0.995677
SNP   | 3306945 | 20551 | 6042 | 0.993824 | 0.998177  | 0.995996

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

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 1022    | 29   | 19   | 0.972407 | 0.982075  | 0.977217
SNP   | 24964   | 315  | 167  | 0.987539 | 0.993355  | 0.990439


## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 115m
call_variants                    | 206m (175m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 71m
total                            | 392m = 6.5 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

(The input BAM is phased already and DeepVariant was run with
`--use_hp_information=true`.)

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel |  501517 | 2984 | 2837 | 0.994085 | 0.994597  | 0.994341
SNP   | 3323836 | 3659 | 2146 | 0.998900 | 0.999355  | 0.999128

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 153m
call_variants                    | 252m (201m with OpenVINO<sup>[*](#vfootnote1)</sup>)
postprocess_variants (with gVCF) | 62m
total                            | 467 m = 7.8 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 504501  | 1330 | 2144 | 0.997364 | 0.995959  | 0.996661
SNP   | 3327495 | 3841 | 1904 | 0.998846 | 0.999428  | 0.999137

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
