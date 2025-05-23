# Runtime and accuracy metrics for Pangenome-aware DeepVariant

## Setup

The runtime and accuracy reported in this page are generated using
`n2-standard-96` GCP instances which has the following configuration:

```bash
GCP instance type: n2-standard-96
CPUs: 96-core (vCPU)
Memory: 384GiB
GPUs: 0
```

## WGS (Illumina)

BAM: We used the VG Giraffe mapped BAM file. The file is available here:
`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
load_gbz_into_shared_memory      | 1m8.16s
make_examples                    | 88m40.61s
call_variants                    | 164m8.08s
postprocess_variants (with gVCF) | 7m19.23s
vcf_stats_report (optional)      | 1m13.14s
total                            | 275m15.26s (4h35m15.26s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502333   | 2168     | 1499     | 0.995703      | 0.997145         | 0.996423        |
| SNP   | 3320003  | 7492     | 5039     | 0.997748      | 0.998485         | 0.998117        |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
load_gbz_into_shared_memory      | 1m8.12s
make_examples                    | 4m57.62s
call_variants                    | 1m4.49s
postprocess_variants (with gVCF) | 0m38.81s
vcf_stats_report (optional)      | 0m5.07s
total                            | 9m22.64s

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1023     | 28       | 16       | 0.973359      | 0.984906         | 0.979098        |
| SNP   | 25006    | 273      | 54       | 0.989201      | 0.997845         | 0.993504        |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](https://github.com/google/deepvariant/blob/r1.9/docs/deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.9/scripts/inference_deepvariant.sh

# WGS-PANGENOME
bash inference_deepvariant.sh --model_preset WGS_PANGENOME

# WGS-PANGENOME
bash inference_deepvariant.sh --model_preset WES_PANGENOME
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant.

The accuracy metrics came from the hap.py program.


