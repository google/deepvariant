# Runtime and accuracy metrics for Pangenome-aware DeepVariant

## How to reproduce the metrics on this page

We report runtime with a n1-highmem-96 machine, and we ran with
`--num_shards 50` to avoid OOM.

This is NOT the fastest or cheapest configuration. We plan to improve the high
memory usage in the future.

## WGS (Illumina)

BAM: We used the VG Giraffe mapped BAM file. You can find it in:
`gs://deepvariant/vg-case-study/HG003.novaseq.pcr-free.35x.vg-1.55.0.bam`

### Runtime

Runtime is on HG003 (all chromosomes).
Reported runtime is an average of 5 runs.

Stage                            | Time (minutes)
-------------------------------- | ------------------
make_examples                    | 170m2.87s
call_variants                    | 522m44.28s
postprocess_variants (with gVCF) | 10m20.44s
vcf_stats_report (optional)      | 9m4.78s
total                            | 721m14.42s (12h1m14.42s)

### Accuracy

hap.py results on HG003 (all chromosomes, using NIST v4.2.1 truth), which was
held out while training.

TODO: Update this table after retraining with pangenome=GBZ. Note that
these numbers are not a correct base for future comparison.

| Type  | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| ----- | -------- | -------- | -------- | ------------- | ---------------- | --------------- |
| INDEL | 501918   | 2583     | 1757     | 0.99488       | 0.996652         | 0.995765        |
| SNP   | 3318596  | 8900     | 8158     | 0.997325      | 0.997549         | 0.997437        |

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

