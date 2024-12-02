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
make_examples                    | 85m58.66s
call_variants                    | 313m49.80s
postprocess_variants (with gVCF) | 7m52.00s
vcf_stats_report (optional)      | 5m48.88s
total                            | 423m6.93s (7h3m6.93s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502338   | 2163     | 1511     | 0.995713      | 0.997122         | 0.996417        |
| SNP   | 3320044  | 7452     | 4735     | 0.99776       | 0.998577         | 0.998168        |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 5m4.22s
call_variants                    | 1m50.67s
postprocess_variants (with gVCF) | 0m38.74s
vcf_stats_report (optional)      | 0m4.91s
total                            | 10m20.44s

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1020     | 31       | 15       | 0.970504      | 0.985782         | 0.978083        |
| SNP   | 25006    | 273      | 54       | 0.989201      | 0.997845         | 0.993504        |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](https://github.com/google/deepvariant/blob/r1.8/docs/deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.8/scripts/inference_deepvariant.sh

# WGS-PANGENOME
bash inference_deepvariant.sh --model_preset WGS_PANGENOME

# WGS-PANGENOME
bash inference_deepvariant.sh --model_preset WES_PANGENOME
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant.

The accuracy metrics came from the hap.py program.


