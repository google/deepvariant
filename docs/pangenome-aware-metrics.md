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

Details of metrics can be found here:

* [Sample sheet](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/pangenome-aware-deepvariant/pangenome_aware_deepvariant_case_study_summaries/sample_sheet.tsv)
* [Multi QC report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/pangenome-aware-deepvariant/pangenome_aware_deepvariant_case_study_summaries/multiqc_report.html)
* [Runtime summary report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/pangenome-aware-deepvariant/pangenome_aware_deepvariant_case_study_summaries/runtimes.md)
* [Accuracy summary report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/pangenome-aware-deepvariant/pangenome_aware_deepvariant_case_study_summaries/happy.summary.md)

Sample sheet contains details of the input files used to generate this report.

Note: Each model type uses different coverages.

## Accuracy

| Model type    | sample   | Type   |   TRUTH.TOTAL |   TRUTH.TP |   TRUTH.FN |   QUERY.TOTAL |   QUERY.FP |   Recall |   Precision |   F1_Score |
|:--------------|:---------|:-------|--------------:|-----------:|-----------:|--------------:|-----------:|---------:|------------:|-----------:|
| wgs-pangenome | HG003    | INDEL  |        504501 |     502667 |       1834 |        940230 |       1263 | 0.996365 |    0.997587 |   0.996976 |
| wgs-pangenome | HG003    | SNP    |       3327496 |    3320057 |       7439 |       4074366 |       4945 | 0.997764 |    0.998514 |   0.998139 |
| wes-pangenome | HG003    | INDEL  |          1051 |       1025 |         26 |          1497 |         14 | 0.975262 |    0.986805 |   0.980999 |
| wes-pangenome | HG003    | SNP    |         25279 |      25007 |        272 |         27679 |         53 | 0.98924  |    0.997885 |   0.993544 |

## Runtime

### Total runtime

| Model type    | sample   | stage                | mean runtime    |
|:--------------|:---------|:---------------------|:----------------|
| wgs-pangenome | HG003    | total                | 4h 8m 47s       |
| wes-pangenome | HG003    | total                | 6m 32s          |

### Runtime by stage

| Model type    | sample   | stage                | mean runtime    |
|:--------------|:---------|:---------------------|:----------------|
| wgs-pangenome | HG003    | make_examples        | 1h 31m 38s      |
| wgs-pangenome | HG003    | call_variants        | 2h 30m 56s      |
| wgs-pangenome | HG003    | postprocess_variants | 6m 12s          |
| wgs-pangenome | HG003    | vcf_stats            | 5m 35s          |
| wgs-pangenome | HG003    | total                | 4h 8m 47s       |
| wes-pangenome | HG003    | make_examples        | 4m 57s          |
| wes-pangenome | HG003    | call_variants        | 1m 4s           |
| wes-pangenome | HG003    | postprocess_variants | 31s             |
| wes-pangenome | HG003    | vcf_stats            | 6s              |
| wes-pangenome | HG003    | total                | 6m 32s          |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
# Get the script.
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.10/scripts/inference_deepvariant.sh

# WGS-PANGENOME
bash inference_deepvariant.sh --model_preset WGS_PANGENOME

# WES-PANGENOME
bash inference_deepvariant.sh --model_preset WES_PANGENOME
```

Runtime metrics are taken from the resulting log after each stage of
DeepVariant.

The accuracy metrics came from the hap.py program.

