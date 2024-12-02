# Runtime and accuracy metrics for all release models

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

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 54m58.62s
call_variants                    | 38m45.29s
postprocess_variants (with gVCF) | 8m22.88s
vcf_stats_report (optional)      | 5m37.52s (optional)
total                            | 113m11.70s (1h53m11.70s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501653   | 2848     | 1289     | 0.994355      | 0.997541         | 0.995945        |
| SNP   | 3306740  | 20756    | 4386     | 0.993762      | 0.998676         | 0.996213        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.8.0/WGS/deepvariant.output.visual_report.html)

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 3m17.64s
call_variants                    | 0m56.36s
postprocess_variants (with gVCF) | 0m39.27s
vcf_stats_report (optional)      | 0m4.93s (optional)
total                            | 5m26.00s

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1020     | 31       | 7        | 0.970504      | 0.993327         | 0.981783        |
| SNP   | 24984    | 295      | 60       | 0.98833       | 0.997604         | 0.992946        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.8.0/WES/deepvariant.output.visual_report.html)

## PacBio (HiFi)

### Updated dataset in release 1.8.0

In release 1.8.0, we have updated the PacBio test data from HG003 Sequel-II to
latest Revio with SPRQ chemistry data to showcase performance on the updated
platform and chemistry. The numbers reported here are generated using the bam
that can be found in:

```bash
gs://deepvariant/pacbio-case-study-testdata/HG003.SPRQ.pacbio.GRCh38.nov2024.bam
```

Which is also available through [here](https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/GIAB_trio/HG003/analysis/GRCh38.m84039_241002_000337_s3.hifi_reads.bc2020.bam).

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | 31m51.00s
call_variants                    | 34m49.62s
postprocess_variants (with gVCF) | 5m28.59s
vcf_stats_report (optional)      | 5m36.49s (optional)
total                            | 86m50.09s (1h26m50.09s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 500955   | 3546     | 3373     | 0.992971      | 0.993555         | 0.993263        |
| SNP   | 3321825  | 5670     | 4263     | 0.998296      | 0.99872          | 0.998508        |


[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.8.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | 53m25.60s
call_variants                    | 55m24.86s
postprocess_variants (with gVCF) | 7m17.83s
vcf_stats_report (optional)      | 6m30.29s (optional)
total                            | 127m56.44s (2h7m56.44s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 452010   | 52491    | 40289    | 0.895955      | 0.920501         | 0.908062        |
| SNP   | 3321452  | 6032     | 3942     | 0.998187      | 0.998815         | 0.998501        |


[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.8.0/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 71m52.43s
call_variants                    | 51m42.37s
postprocess_variants (with gVCF) | 4m6.13s
vcf_stats_report (optional)      | 5m18.39s (optional)
total                            | 151m34.49s (2h31m34.49s)

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503109   | 1392     | 2636     | 0.997241      | 0.995022         | 0.99613         |
| SNP   | 3324179  | 3316     | 2049     | 0.999003      | 0.999384         | 0.999194        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.8.0/HYBRID/deepvariant.output.visual_report.html)

## Inspect outputs that produced the metrics above

The DeepVariant VCFs, gVCFs, and hap.py evaluation outputs are available at:

```
gs://deepvariant/case-study-outputs
```

You can also inspect them in a web browser here:
https://42basepairs.com/browse/gs/deepvariant/case-study-outputs

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.8/scripts/inference_deepvariant.sh

# WGS
bash inference_deepvariant.sh --model_preset WGS

# WES
bash inference_deepvariant.sh --model_preset WES

# PacBio
bash inference_deepvariant.sh --model_preset PACBIO

# ONT_R104
bash inference_deepvariant.sh --model_preset ONT_R104

# Hybrid
bash inference_deepvariant.sh --model_preset HYBRID_PACBIO_ILLUMINA
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant. The runtime numbers reported above are the average of 5 runs each.
The accuracy metrics come from the hap.py summary.csv output file.
The runs are deterministic so all 5 runs produced the same output.

[CPU instance with 96 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform

