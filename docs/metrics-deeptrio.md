# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -----------------
make_examples                    | ~460m
call_variants for HG002-4        | ~475m
call_variants for HG003          | ~475m
call_variants for HG004          | ~475m
postprocess_variants (parallel)  | ~110m
total                            | ~1995m = 33.3 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

#### HG002:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 523006  | 2463  | 839  | 0.995313 | 0.998464  | 0.996886
SNP   | 3347854 | 17273 | 3057 | 0.994867 | 0.999088  | 0.996973

#### HG003:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 501697  | 2804  | 1234 | 0.994442 | 0.997646  | 0.996041
SNP   | 3308641 | 18855 | 3373 | 0.994334 | 0.998982  | 0.996652

#### HG004:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 507681  | 2838  | 1182 | 0.994441 | 0.997774  | 0.996104
SNP   | 3326961 | 19649 | 3232 | 0.994129 | 0.999030  | 0.996573

## PacBio (HiFi)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | -------------------
make_examples                    | ~840m
call_variants for HG002          | ~400m
call_variants for HG003          | ~400m
call_variants for HG004          | ~400m
postprocess_variants (parallel)  | ~80m
total                            | ~2120m = 35.3 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth). HG002/HG003/HG004 trio is used for training a PacBio model, only chr20,
chr21, and chr22 are held out. Accuracy metrics are reported for the whole
genome for consistency.

#### HG002:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 524058  | 1411 | 1865 | 0.997315 | 0.99661   | 0.996962
SNP   | 3361591 | 3536 | 980  | 0.998949 | 0.999709  | 0.999329

#### HG003:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 502477  | 2024 | 2333 | 0.995988 | 0.995578  | 0.995783
SNP   | 3323668 | 3827 | 1509 | 0.998850 | 0.999547  | 0.999198

#### HG004:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 508346  | 2173 | 2471 | 0.995744 | 0.995375  | 0.995559
SNP   | 3342956 | 3654 | 1351 | 0.998908 | 0.999596  | 0.999252

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Wall time (minutes)
-------------------------------- | --------------
make_examples                    | ~18m
call_variants for HG002          | ~7m
call_variants for HG003          | ~7m
call_variants for HG004          | ~7m
postprocess_variants (parallel)  | ~2m
total                            | ~41m

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2.1
truth).

#### HG002:

Type  | # TP  | # FN | # FP | Recall  | Precision | F1_Score
----- | ----- | ---- | ---- | ------- | --------- | --------
Indel | 1059  | 31   | 12   | 0.97156 | 0.989011  | 0.980208
SNP   | 25165 | 249  | 31   | 0.990202| 0.998770  | 0.994468

#### HG003:

Type  | # TP  | # FN | # FP | Recall   | Precision | F1_Score
----- | ----- | ---- | ---- | -------- | --------- | --------
Indel | 1022  | 29   | 11   | 0.972407 | 0.989564  | 0.980910
SNP   | 25038 | 241  | 28   | 0.990466 | 0.998883  | 0.994657


#### HG004:

Type  | # TP  | # FN | # FP | Recall   | Precision | F1_Score
----- | ----- | ---- | ---- | -------- | --------- | --------
Indel | 1057  | 25   | 20   | 0.976895 | 0.981982  | 0.979432
SNP   | 24945 | 232  | 17   | 0.990785 | 0.999319  | 0.995034


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
