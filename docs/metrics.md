# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 94m
call_variants                    | 215m
postprocess_variants (with gVCF) | 86m
total                            | 395m = 6.6 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2 truth), which was held
out while training.

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel |  501841 |  3069 | 1389 | 0.993922 | 0.997351  | 0.995634
SNP   | 3310730 | 20760 | 6202 | 0.993769 | 0.998131  | 0.995945

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

hap.py results on HG003 (all chromosomes, using NIST v4.2 truth), which was held
out while training.

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 1025    | 28   | 19   | 0.973409 | 0.982126  | 0.977748
SNP   | 25005   | 319  | 169  | 0.987403 | 0.993287  | 0.990337


## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 115m
call_variants                    | 195m
postprocess_variants (with gVCF) | 72m
total                            | 382m = 6.4 hours

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2 truth), which was held
out while training.

(The input BAM is phased already and DeepVariant
was run with `--sort_by_haplotypes=true --parse_sam_aux_fields=true`.)

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel |  501468 | 3442 | 3461 | 0.993183 | 0.993416  | 0.993300
SNP   | 3327592 | 3898 | 2535 | 0.998830 | 0.999239  | 0.999035

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 145m
call_variants                    | 230m
postprocess_variants (with gVCF) | 58m
total                            | 433 m = 7.2 hours

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2 truth), which was held out
while training the hybrid model.

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 503570  | 1340 | 2149 | 0.997346 | 0.995953  | 0.996649
SNP   | 3327590 | 3900 | 1934 | 0.998829 | 0.999419  | 0.999124

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration. For more scalable execution
of DeepVariant see the [External Solutions] section.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# WGS (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wgs_case_study_docker.sh
bash inference_wgs.sh

# WES (should take about 20 minutes)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wes_case_study_docker.sh
bash run_wes_case_study_docker.sh

# PacBio (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.1/scripts/inference_pacbio.sh
bash inference_pacbio.sh

# Hybrid (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_hybrid_pacbio_illumina_case_study_docker.sh
bash run_hybrid_pacbio_illumina_case_study_docker.sh
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant, and the accuracy metrics come from the hap.py summary.csv output
file.

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
