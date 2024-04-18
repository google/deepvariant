# Runtime and accuracy metrics for Pangenome-aware DeepVariant

## WGS (Illumina)

BAM: We used the VG Giraffe mapped BAM file. You can find it in:
`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 265m36.96s
call_variants                    | 504m29.99s
postprocess_variants (with gVCF) | 36m40.23s
vcf_stats_report (optional)      | 8m1.51s
total                            | 829m58.70s (13h49m58.70s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 502405   | 2096     | 1427     | 0.995845      | 0.997283         | 0.996564        |
| SNP   | 3319209  | 8287     | 3032     | 0.99751       | 0.999088         | 0.998298        |

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

