# Runtime and accuracy metrics for all release models

## WGS (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 90m39.39s
call_variants                    | 185m21.95s
postprocess_variants (with gVCF) | 22m30.41s
vcf_stats_report (optional)      | 7m20.41s
total                            | 313m52.18s (5h13m52.18s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501676   | 2825     | 1267     | 0.9944        | 0.997583         | 0.995989        |
| SNP   | 3306782  | 20714    | 4284     | 0.993775      | 0.998707         | 0.996235        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.7.0/WGS/deepvariant.output.visual_report.html)

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | 5m52.94s
call_variants                    | 1m18.60s
postprocess_variants (with gVCF) | 0m42.01s
vcf_stats_report (optional)      | 0m6.57s
total                            | 8m30.47s

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 1020     | 31       | 12       | 0.970504      | 0.988615         | 0.979476        |
| SNP   | 24982    | 297      | 64       | 0.988251      | 0.997445         | 0.992827        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.7.0/WES/deepvariant.output.visual_report.html)

## PacBio (HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | 87m9.39s
call_variants                    | 212m12.40s
postprocess_variants (with gVCF) | 25m1.08s
vcf_stats_report (optional)      | 8m12.02s
total                            | 343m28.28s (5h43m28.28s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

Starting from v1.4.0, users don't need to phase the BAMs first, and only need
to run DeepVariant once.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501504   | 2997     | 2760     | 0.994059      | 0.994744         | 0.994402        |
| SNP   | 3324308  | 3187     | 1508     | 0.999042      | 0.999547         | 0.999295        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.7.0/PACBIO/deepvariant.output.visual_report.html)

## ONT_R104

### Runtime

Runtime is on HG003 reads (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | --------------------
make_examples                    | 167m43.26s
call_variants                    | 264m34.06s
postprocess_variants (with gVCF) | 28m21.94s
vcf_stats_report (optional)      | 8m44.50s
total                            | 472m31.96s (7h52m31.96s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1
truth), which was held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 441622   | 62879    | 41272    | 0.875364      | 0.917464         | 0.89592         |
| SNP   | 3314062  | 13433    | 8137     | 0.995963      | 0.997552         | 0.996757        |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.7.0/ONT_R104/deepvariant.output.visual_report.html)

## Hybrid (Illumina + PacBio HiFi)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -------------------
make_examples                    | 123m45.96s
call_variants                    | 188m42.36s
postprocess_variants (with gVCF) | 18m32.12s
vcf_stats_report (optional)      | 7m14.30s
total                            | 350m25.79s (5h50m25.79s)

### Accuracy

Evaluating on HG003 (all chromosomes, using NIST v4.2.1 truth), which was held
out while training the hybrid model.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 503013   | 1488     | 2742     | 0.997051      | 0.994828         | 0.995938        |
| SNP   | 3323633  | 3862     | 2259     | 0.998839      | 0.999321         | 0.99908         |

[See VCF stats report.](https://storage.googleapis.com/deepvariant/visual_reports/DeepVariant/1.7.0/HYBRID/deepvariant.output.visual_report.html)

## Inspect outputs that produced the metrics above

The DeepVariant VCFs, gVCFs, and hap.py evaluation outputs are available at:

```
gs://deepvariant/case-study-outputs
```

You can also inspect them in a web browser here:
https://42basepairs.com/browse/gs/deepvariant/case-study-outputs

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.7/scripts/inference_deepvariant.sh

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

[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform

