# Runtime and accuracy metrics for Pangenome-aware DeepVariant

## WGS (Illumina)

BAM: We used the VG Giraffe mapped BAM file. You can find it in:
`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 330m38.83s
call_variants                    | 534m22.24s
postprocess_variants (with gVCF) | 40m42.87s
vcf_stats_report (optional)      | 8m47.20s
total                            | 929m51.74s (15h29m51.74s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502505   | 1996     | 1367     | 0.996044      | 0.997397         | 0.99672         |
| SNP   | 3319187  | 8309     | 4395     | 0.997503      | 0.998678         | 0.99809         |

## WES (Illumina)

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | -----------------
make_examples                    | TODO
call_variants                    | TODO
postprocess_variants (with gVCF) | TODO
vcf_stats_report (optional)      | TODO
total                            | TODO

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

TODO

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 64 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

[CPU instance with 64 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform

