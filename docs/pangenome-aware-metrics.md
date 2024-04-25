# Runtime and accuracy metrics for Pangenome-aware DeepVariant

## WGS (Illumina)

BAM: We used the VG Giraffe mapped BAM file. You can find it in:
`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 269m21.26s
call_variants                    | 525m21.01s
postprocess_variants (with gVCF) | 39m49.49s
vcf_stats_report (optional)      | 8m43.36s
total                            | 858m25.33s (14h18m25.33s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502516   | 1985     | 1315     | 0.996065      | 0.997496         | 0.99678         |
| SNP   | 3320314  | 7182     | 3580     | 0.997842      | 0.998924         | 0.998382        |

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

