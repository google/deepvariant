# Runtime and accuracy metrics for all release models

## WGS (Illumina)

## WES (Illumina)

## PacBio (HiFi)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 138m
call_variants                    | 265m
postprocess_variants (with gVCF) | 59m
total                            | 462 m = 7.1 hours

### Accuracy

Evaluating on HG003, which was held out while training the hybrid model.

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 503570  | 1340 | 2149 | 0.997346 | 0.995953  | 0.996649
SNP   | 3327590 | 3900 | 1934 | 0.998829 | 0.999419  | 0.999124

## How to reproduce the metrics on this page

Create an instance with 64 CPUs:

```
gcloud compute instances create deepvariant-hybrid-case-study \
--scopes "compute-rw,storage-full,cloud-platform" \
--image-family "ubuntu-1604-lts" \
--image-project "ubuntu-os-cloud" \
--zone "us-west1-b" \
--machine-type "custom-64-131072" \
--min-cpu-platform "Intel Skylake" \
--boot-disk-size "300" \
--boot-disk-type "pd-ssd"
```

Log in to the newly created instance:

```
gcloud compute ssh "deepvariant-hybrid-case-study" --zone "us-west1-b"
```

Download and run any of the following case study scripts:

```
# WGS (should take about 6 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wgs_case_study_docker.sh
bash run_wgs_case_study_docker.sh

# WES (should take about 20 minutes)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_wes_case_study_docker.sh
bash run_wes_case_study_docker.sh

# PacBio (should take about 6 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_pacbio_case_study_docker.sh
bash run_pacbio_case_study_docker.sh

# Hybrid (should take about 7 hours)
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.0/scripts/run_hybrid_pacbio_illumina_case_study_docker.sh
bash run_hybrid_pacbio_illumina_case_study_docker.sh
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant, and the accuracy metrics come from the hap.py summary.csv output
file.
