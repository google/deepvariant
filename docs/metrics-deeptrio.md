# DeepTrio runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | ~440m
call_variants for HG002          | ~500m
call_variants for HG003          | ~500m
call_variants for HG004          | ~500m
postprocess_variants (with gVCF) | ~110m
total                            | ~2050m ~ 34 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2
truth), which was held out while training.

#### HG002:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 523008  | 2458  | 840  | 0.995322 | 0.998462  | 0.996890
SNP   | 3348025 | 17315 | 3054 | 0.994855 | 0.999089  | 0.996968

#### HG003:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 502068  | 2842  | 1236 | 0.994371 | 0.997644  | 0.996005
SNP   | 3312429 | 19061 | 3449 | 0.994279 | 0.998960  | 0.996614

#### HG004:

Type  | # TP    | # FN  | # FP | Recall   | Precision | F1_Score
----- | ------- | ----- | ---- | -------- | --------- | --------
Indel | 508628  | 2898  | 1198 | 0.994335 | 0.997748  | 0.996038
SNP   | 3335516 | 20082 | 3340 | 0.994015 | 0.999000  | 0.996502

## PacBio (HiFi)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | ~840m
call_variants for HG002          | ~400m
call_variants for HG003          | ~400m
call_variants for HG004          | ~400m
postprocess_variants (with gVCF) | ~84m
total                            | ~2124m = 35.4 hours

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2
truth). HG002/HG003/HG004 trio is used for training a PacBio model, only chr20,
chr21, and chr22 are held out. Accuracy metrics are reported for the whole
genome for consistency.

#### HG002:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 524061  | 1405 | 1864 | 0.997326 | 0.996612  | 0.996969
SNP   | 3361815 | 3525 | 969  | 0.998953 | 0.999712  | 0.999332

#### HG003:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 502884  | 2026 | 2330 | 0.995987 | 0.995587  | 0.995787
SNP   | 3327628 | 3862 | 1522 | 0.998841 | 0.999543  | 0.999192

#### HG004:

Type  | # TP    | # FN | # FP | Recall   | Precision | F1_Score
----- | ------- | ---- | ---- | -------- | --------- | --------
Indel | 509342  | 2184 | 2483 | 0.995730 | 0.995361  | 0.995546
SNP   | 3351890 | 3708 | 1389 | 0.998895 | 0.999586  | 0.999240

## Whole Exome Sequencing (Illumina)

### Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).

Stage                            | Time (minutes)
-------------------------------- | --------------
make_examples                    | ~18m
call_variants for HG002          | ~6m
call_variants for HG003          | ~6m
call_variants for HG004          | ~6m
postprocess_variants (with gVCF) | ~2m
total                            | ~38m

### Accuracy

hap.py results on HG002/HG003/HG004 trio (all chromosomes, using NIST v4.2
truth).

#### HG002:

Type  | # TP  | # FN | # FP | Recall  | Precision | F1_Score
----- | ----- | ---- | ---- | ------- | --------- | --------
Indel | 1059  | 31   | 12   | 0.97156 | 0.989011  | 0.980208
SNP   | 25179 | 251  | 31   | 0.99013 | 0.998770  | 0.994431

#### HG003:

Type  | # TP  | # FN | # FP | Recall   | Precision | F1_Score
----- | ----- | ---- | ---- | -------- | --------- | --------
Indel | 1025  | 28   | 11   | 0.973409 | 0.989593  | 0.981435
SNP   | 25079 | 245  | 28   | 0.990325 | 0.998885  | 0.994587

#### HG004:

Type  | # TP  | # FN | # FP | Recall   | Precision | F1_Score
----- | ----- | ---- | ---- | -------- | --------- | --------
Indel | 1061  | 25   | 20   | 0.976980 | 0.982047  | 0.979507
SNP   | 25008 | 241  | 17   | 0.990455 | 0.999321  | 0.994868

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
