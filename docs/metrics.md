# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 85m
call_variants                    | 231m
postprocess_variants (with gVCF) | 87m
total                            | 431m = 6.7 hours

### Accuracy

hap.py results on HG002 (all chromosomes), using NIST v4.2 truth.

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 522259  | 3207  | 1187 | 0.993897 | 0.997825  | 0.995857
SNP   | 3345988 | 19352 | 3955 | 0.994250 | 0.998820  | 0.996530


## WES (Illumina)

### Runtime

Runtime is on HG002 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 18m
call_variants                    | 2m
postprocess_variants (with gVCF) | 1m
total                            | 21m

### Accuracy

hap.py results on HG002 (all chromosomes), using NIST v4.1 truth.

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 2896    | 124  | 81   | 0.958940 | 0.973134  | 0.965985
SNP   | 38180   | 396  | 116  | 0.989735 | 0.996975  | 0.993341


## PacBio (HiFi)

### Runtime

Runtime is on HG002 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 152m
call_variants                    | 204m
postprocess_variants (with gVCF) | 65m
total                            | 421 = 7 hours

### Accuracy

hap.py results on HG002 (chr20), using NIST v4.2 truth.

(The input BAM is haplotagged already and DeepVariant
was run with `--sort_by_haplotypes=true --parse_sam_aux_fields=true`.)

Type  | # TP  | # FN | # FP | Recall   | Precision | F1_Score
----- | ----- | ---- | ---- | -------- | --------- | --------
Indel | 11142 | 114  | 111  | 0.989872 | 0.990483  | 0.990177
SNP   | 71273 | 60   | 15   | 0.999159 | 0.999790  | 0.999474

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 136m
call_variants                    | 259m
postprocess_variants (with gVCF) | 60m
total                            | 455 m = 7.6 hours

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
bash run_wgs_case_study_docker.sh

# WES (should take about 20 minutes)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wes_case_study_docker.sh
bash run_wes_case_study_docker.sh

# PacBio (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_pacbio_case_study_docker.sh
bash run_pacbio_case_study_docker.sh

# Hybrid (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_hybrid_pacbio_illumina_case_study_docker.sh
bash run_hybrid_pacbio_illumina_case_study_docker.sh
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant, and the accuracy metrics come from the hap.py summary.csv output
file.

[External Solutions]: https://github.com/google/deepvariant#external-solutions
[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
