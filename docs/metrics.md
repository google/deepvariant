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

Details of metrics can be found here:

* [Sample sheet](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deepvariant_case_study_summaries/sample_sheet.tsv)
* [Multi QC report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deepvariant_case_study_summaries/multiqc_report.html)
* [Runtime summary report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deepvariant_case_study_summaries/runtimes.md)
* [Accuracy summary report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deepvariant_case_study_summaries/happy.summary.md)

Sample sheet contains details of the input files used to generate this report.

Note: Each model type uses different coverages.

## Accuracy

Below we report full genome accuracy as reported using
[hap.py](https://github.com/Illumina/hap.py) for our models.

| Model type             | sample   | Type   |   TRUTH.TOTAL |   TRUTH.TP |   TRUTH.FN |   QUERY.TOTAL |   QUERY.FP |   Recall |   Precision |   F1_Score |
|:-----------------------|:---------|:-------|--------------:|-----------:|-----------:|--------------:|-----------:|---------:|------------:|-----------:|
| wgs                    | HG003    | INDEL  |        504501 |     501594 |       2907 |        937937 |       1190 | 0.994238 |    0.997729 |   0.99598  |
| wgs                    | HG003    | SNP    |       3327496 |    3306720 |      20776 |       3817962 |       4880 | 0.993756 |    0.998527 |   0.996136 |
| wes                    | HG003    | INDEL  |          1051 |       1024 |         27 |          1485 |          8 | 0.97431  |    0.992417 |   0.98328  |
| wes                    | HG003    | SNP    |         25279 |      24983 |        296 |         27709 |         60 | 0.988291 |    0.997604 |   0.992926 |
| pacbio                 | HG003    | INDEL  |        504501 |     501567 |       2934 |        989958 |       3057 | 0.994184 |    0.994162 |   0.994173 |
| pacbio                 | HG003    | SNP    |       3327495 |    3321765 |       5730 |       4329942 |       4125 | 0.998278 |    0.998761 |   0.99852  |
| ont-r104               | HG003    | INDEL  |        504501 |     460355 |      44146 |        830072 |      25676 | 0.912496 |    0.948695 |   0.930243 |
| ont-r104               | HG003    | SNP    |       3327495 |    3321799 |       5696 |       4400475 |       4611 | 0.998288 |    0.998615 |   0.998451 |
| rnaseq                 | HG005    | INDEL  |           188 |        151 |         37 |           285 |         36 | 0.803191 |    0.810526 |   0.806842 |
| rnaseq                 | HG005    | SNP    |         11349 |      10656 |        693 |         12336 |        391 | 0.938937 |    0.964599 |   0.951595 |
| hybrid-pacbio-illumina | HG003    | INDEL  |        504501 |     503264 |       1237 |        998274 |       2052 | 0.997548 |    0.996129 |   0.996838 |
| hybrid-pacbio-illumina | HG003    | SNP    |       3327495 |    3324021 |       3474 |       4068058 |       1856 | 0.998956 |    0.999442 |   0.999199 |

## Runtime

Each case study was run 5x times and the runtimes were averaged. Here we report
the mean runtime in seconds, the standard deviation of runtimes, and a duration
format (`mean_runtime`; hours, minutes, seconds).

### Total runtime only

| Model type             | sample   | stage                | mean runtime    | total_runs |
|:-----------------------|:---------|:---------------------|:----------------|-----------:|
| wgs                    | HG003    | total                | 1h 8m 58s       |          5 |
| exome                    | HG003    | total                | 4m 11s          |          5 |
| pacbio                 | HG003    | total                | 1h 2m 17s       |          5 |
| ont-r104               | HG003    | total                | 1h 43m 18s      |          5 |
| rnaseq                 | HG005    | total                | 9m 1s           |          5 |
| hybrid-pacbio-illumina | HG003    | total                | 2h 14m 54s      |          5 |

### Detailed runtime

| Model type             | sample   | stage                | mean runtime    | total_runs |
|:-----------------------|:---------|:---------------------|:----------------|-----------:|
| wgs                    | HG003    | make_examples        | 46m 15s         |          5 |
| wgs                    | HG003    | call_variants        | 15m 58s         |          5 |
| wgs                    | HG003    | postprocess_variants | 6m 45s          |          5 |
| wgs                    | HG003    | vcf_stats            | 5m 17s          |          5 |
| wgs                    | HG003    | total                | 1h 8m 58s       |          5 |
| exome                  | HG003    | make_examples        | 3m 6s           |          5 |
| exome                  | HG003    | call_variants        | 34s             |          5 |
| exome                  | HG003    | postprocess_variants | 30s             |          5 |
| exome                  | HG003    | vcf_stats            | 6s              |          5 |
| exome                  | HG003    | total                | 4m 11s          |          5 |
| pacbio                 | HG003    | make_examples        | 37m 4s          |          5 |
| pacbio                 | HG003    | call_variants        | 18m 28s         |          5 |
| pacbio                 | HG003    | postprocess_variants | 6m 45s          |          5 |
| pacbio                 | HG003    | vcf_stats            | 5m 46s          |          5 |
| pacbio                 | HG003    | total                | 1h 2m 17s       |          5 |
| ont-r104               | HG003    | make_examples        | 56m 4s          |          5 |
| ont-r104               | HG003    | call_variants        | 32m 52s         |          5 |
| ont-r104               | HG003    | postprocess_variants | 14m 21s         |          5 |
| ont-r104               | HG003    | vcf_stats            | 7m 23s          |          5 |
| ont-r104               | HG003    | total                | 1h 43m 18s      |          5 |
| rnaseq                 | HG005    | make_examples        | 7m 31s          |          5 |
| rnaseq                 | HG005    | call_variants        | 25s             |          5 |
| rnaseq                 | HG005    | postprocess_variants | 1m 4s           |          5 |
| rnaseq                 | HG005    | vcf_stats            | 5s              |          5 |
| rnaseq                 | HG005    | total                | 9m 1s           |          5 |
| hybrid-pacbio-illumina | HG003    | make_examples        | 1h 54s          |          5 |
| hybrid-pacbio-illumina | HG003    | call_variants        | 1h 10m 4s       |          5 |
| hybrid-pacbio-illumina | HG003    | postprocess_variants | 3m 55s          |          5 |
| hybrid-pacbio-illumina | HG003    | vcf_stats            | 5m 3s           |          5 |
| hybrid-pacbio-illumina | HG003    | total                | 2h 14m 54s      |          5 |

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
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.10/scripts/inference_deepvariant.sh

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
The accuracy metrics come from the hap.py summary.csv output file. The runs are
deterministic so all 5 runs produced the same output.

