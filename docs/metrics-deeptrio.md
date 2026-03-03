# DeepTrio runtime and accuracy metrics for all release models

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

* [Sample sheet](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deeptrio/deeptrio_case_study_summaries/sample_sheet.tsv)
* [Multi QC report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deeptrio/deeptrio_case_study_summaries/multiqc_report.html)
* [Accuracy summary report](https://storage.googleapis.com/deepvariant/case-study-outputs/1.10.0/deeptrio/deeptrio_case_study_summaries/happy.summary.md)

Sample sheet contains details of the input files used to generate this report.

Note: Each model type uses different coverages.

## Accuracy

We report accuracy on `chr20` only.

| Model type     | sample   | Type   |   TRUTH.TOTAL |   TRUTH.TP |   TRUTH.FN |   QUERY.TOTAL |   QUERY.FP |   Recall |   Precision |   F1_Score |
|:---------------|:---------|:-------|--------------:|-----------:|-----------:|--------------:|-----------:|---------:|------------:|-----------:|
| wgs-chr20      | HG002    | INDEL  |         11256 |      11204 |         52 |         21103 |         17 | 0.99538  |    0.998541 |   0.996958 |
| wgs-chr20      | HG002    | SNP    |         71333 |      71087 |        246 |         88094 |         39 | 0.996551 |    0.999452 |   0.998    |
| wgs-chr20      | HG003    | INDEL  |         10628 |      10571 |         57 |         20941 |         28 | 0.994637 |    0.997467 |   0.99605  |
| wgs-chr20      | HG003    | SNP    |         70166 |      69975 |        191 |         85456 |         55 | 0.997278 |    0.999215 |   0.998246 |
| wgs-chr20      | HG004    | INDEL  |         11000 |      10945 |         55 |         21348 |         29 | 0.995    |    0.997469 |   0.996233 |
| wgs-chr20      | HG004    | SNP    |         71659 |      71445 |        214 |         86665 |         49 | 0.997014 |    0.999315 |   0.998163 |
| wes-chr20      | HG002    | INDEL  |            34 |         34 |          0 |            41 |          0 | 1        |    1        |   1        |
| wes-chr20      | HG002    | SNP    |           672 |        670 |          2 |           708 |          0 | 0.997024 |    1        |   0.99851  |
| wes-chr20      | HG003    | INDEL  |            29 |         29 |          0 |            46 |          0 | 1        |    1        |   1        |
| wes-chr20      | HG003    | SNP    |           685 |        683 |          2 |           703 |          0 | 0.99708  |    1        |   0.998538 |
| wes-chr20      | HG004    | INDEL  |            33 |         32 |          1 |            45 |          1 | 0.969697 |    0.969697 |   0.969697 |
| wes-chr20      | HG004    | SNP    |           679 |        676 |          3 |           707 |          0 | 0.995582 |    1        |   0.997786 |
| pacbio-chr20   | HG002    | INDEL  |         11256 |      11228 |         28 |         23150 |         57 | 0.997512 |    0.995154 |   0.996332 |
| pacbio-chr20   | HG002    | SNP    |         71333 |      71311 |         22 |        109348 |         19 | 0.999692 |    0.999734 |   0.999713 |
| pacbio-chr20   | HG003    | INDEL  |         10628 |      10580 |         48 |         23643 |         70 | 0.995484 |    0.993715 |   0.994599 |
| pacbio-chr20   | HG003    | SNP    |         70166 |      70147 |         19 |        118404 |         37 | 0.999729 |    0.999473 |   0.999601 |
| pacbio-chr20   | HG004    | INDEL  |         11000 |      10958 |         42 |         23959 |         58 | 0.996182 |    0.994968 |   0.995574 |
| pacbio-chr20   | HG004    | SNP    |         71659 |      71621 |         38 |        118676 |         26 | 0.99947  |    0.999638 |   0.999554 |

## Runtime

Runtime is on HG002/HG003/HG004 (all chromosomes).
Reported runtime is an average of 5 runs.

### Total runtime

| Model type       | mean runtime    |
|:-----------------|:----------------|
| wgs              | 4h 37m 19s      |
| pacbio           | 6h 14m 8s       |
| wes              | 14m 27s         |

### Runtime by stage
| Model type| stage                        | mean runtime  |
|-----------|------------------------------|---------------|
| wgs       | make_examples                | 2h 54m 56s    |
| wgs       | call_variants_child          | 20m 7s        |
| wgs       | call_variants_parent1        | 21m 56s       |
| wgs       | call_variants_parent2        | 21m 42s       |
| wgs       | postprocess_variants_child   | 7m 50s        |
| wgs       | postprocess_variants_parent1 | 8m 15s        |
| wgs       | postprocess_variants_parent2 | 8m 22s        |
| wgs       | vcf_stats_report_child       | 6m 26s        |
| wgs       | vcf_stats_report_parent1     | 6m 30s        |
| wgs       | vcf_stats_report_parent2     | 6m 39s        |
| wes       | make_examples                | 7m 10s        |
| wes       | call_variants_child          | 2m 9s         |
| wes       | call_variants_parent1        | 2m 11s        |
| wes       | call_variants_parent2        | 2m 11s        |
| wes       | postprocess_variants_child   | 9s            |
| wes       | postprocess_variants_parent1 | 9s            |
| wes       | postprocess_variants_parent2 | 10s           |
| wes       | vcf_stats_report_child       | 6s            |
| wes       | vcf_stats_report_parent1     | 6s            |
| wes       | vcf_stats_report_parent2     | 6s            |
| pacbio    | make_examples                | 3h 40m 15s    |
| pacbio    | call_variants_child          | 28m 42s       |
| pacbio    | call_variants_parent1        | 35m 23s       |
| pacbio    | call_variants_parent2        | 35m 35s       |
| pacbio    | postprocess_variants_child   | 7m 33s        |
| pacbio    | postprocess_variants_parent1 | 7m 45s        |
| pacbio    | postprocess_variants_parent2 | 8m 2s         |
| pacbio    | vcf_stats_report_child       | 6m 55s        |
| pacbio    | vcf_stats_report_parent1     | 7m 1s         |
| pacbio    | vcf_stats_report_parent2     | 7m 3s         |

## How to reproduce the metrics on this page

For simplicity and consistency, we report runtime with a
[CPU instance with 96 CPUs](deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform)
For bigger datasets (WGS and PACBIO), we used bigger disk size (900G).
This is NOT the fastest or cheapest configuration.

Use `gcloud compute ssh` to log in to the newly created instance.

Download and run any of the following case study scripts:

```
curl -O https://raw.githubusercontent.com/google/deepvariant/r1.10/scripts/inference_deeptrio.sh

# WGS
bash inference_deeptrio.sh --model_preset WGS

# WES
bash inference_deeptrio.sh --model_preset WES

# PacBio
bash inference_deeptrio.sh --model_preset PACBIO

```

Runtime metrics are taken from the resulting log after each stage of
DeepTrio. The runtime numbers reported above are the average of 5 runs each.
The accuracy metrics come from the hap.py summary.csv output file.
The runs are deterministic so all 5 runs produced the same output.

[CPU instance with 96 CPUs]: deepvariant-details.md#command-for-a-cpu-only-machine-on-google-cloud-platform
