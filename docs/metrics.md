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

Reported values are based on evaluations of HG003.

## Accuracy

Below we report full genome accuracy as reported using
[hap.py](https://github.com/Illumina/hap.py).

model_type             | Type  | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | Recall   | Precision | F1_Score
:--------------------- |:----- | ----------: | -------: | -------: | ----------: | -------: | -------: | --------: | -------:
wgs                    | INDEL | 504501      | 501594   | 2907     | 937937      | 1190     | 0.994238 | 0.997729  | 0.99598
wgs                    | SNP   | 3327496     | 3306720  | 20776    | 3817962     | 4880     | 0.993756 | 0.998527  | 0.996136
exome                  | INDEL | 1051        | 1024     | 27       | 1485        | 8        | 0.97431  | 0.992417  | 0.98328
exome                  | SNP   | 25279       | 24983    | 296      | 27709       | 60       | 0.988291 | 0.997604  | 0.992926
pacbio                 | INDEL | 504501      | 501598   | 2903     | 986955      | 2949     | 0.994246 | 0.994368  | 0.994307
pacbio                 | SNP   | 3327495     | 3321742  | 5753     | 4331772     | 4107     | 0.998271 | 0.998767  | 0.998519
ont-r104               | INDEL | 504501      | 463074   | 41427    | 895345      | 35116    | 0.917885 | 0.931685  | 0.924733
ont-r104               | SNP   | 3327495     | 3321037  | 6458     | 4408429     | 5729     | 0.998059 | 0.998279  | 0.998169
hybrid-pacbio-illumina | INDEL | 504501      | 503264   | 1237     | 998274      | 2052     | 0.997548 | 0.996129  | 0.996838
hybrid-pacbio-illumina | SNP   | 3327495     | 3324021  | 3474     | 4068058     | 1856     | 0.998956 | 0.999442  | 0.999199

## Runtime

Each case study was run 5x times and the runtimes were averaged. Here we report
the mean runtime in seconds, the standard deviation of runtimes, and a duration
format (`mean_hruntime`; hours, minutes, seconds).

model_type             | stage                | mean_runtime (s) | std_runtime | mean_hruntime
:--------------------- | :------------------- | ---------------: | ----------: | :------------
wgs                    | make_examples        | 2887.1           | 68.658      | 48m 7s
wgs                    | call_variants        | 939.88           | 19.599      | 15m 39s
wgs                    | postprocess_variants | 403.37           | 3.327       | 6m 43s
wgs                    | vcf_stats            | 317.07           | 1.123       | 5m 17s
wgs                    | total                | 4230.35          |             | 1h 10m 30s
exome                  | make_examples        | 176.57           | 2.153       | 2m 56s
exome                  | call_variants        | 33.28            | 0.224       | 33s
exome                  | postprocess_variants | 29.28            | 0.465       | 29s
exome                  | vcf_stats            | 4.95             | 0.046       | 4s
exome                  | total                | 239.13           |             | 3m 59s
pacbio                 | make_examples        | 2036.71          | 104.087     | 33m 56s
pacbio                 | call_variants        | 697.31           | 61.092      | 11m 37s
pacbio                 | postprocess_variants | 291.27           | 6.432       | 4m 51s
pacbio                 | vcf_stats            | 340.26           | 11.488      | 5m 40s
pacbio                 | total                | 3025.29          |             | 50m 25s
ont-r104               | make_examples        | 3042.24          | 20.359      | 50m 42s
ont-r104               | call_variants        | 3286.89          | 104.469     | 54m 46s
ont-r104               | postprocess_variants | 669.59           | 5.558       | 11m 9s
ont-r104               | vcf_stats            | 444.71           | 10.684      | 7m 24s
ont-r104               | total                | 6998.72          |             | 1h 56m 38s
hybrid-pacbio-illumina | make_examples        | 3648.28          | 34.422      | 1h 48s
hybrid-pacbio-illumina | call_variants        | 4215.97          | 314.295     | 1h 10m 15s
hybrid-pacbio-illumina | postprocess_variants | 235.97           | 2.797       | 3m 55s
hybrid-pacbio-illumina | vcf_stats            | 305.55           | 1.529       | 5m 5s
hybrid-pacbio-illumina | total                | 8100.22          |             | 2h 15m

**Total Runtime**

The total rows are summarized below as well:

uid                    | sample | mean_hruntime
:--------------------- | :----- | :------------
wgs                    | HG003  | 1h 10m 30s
exome                  | HG003  | 3m 59s
pacbio                 | HG003  | 50m 25s
ont-r104               | HG003  | 1h 56m 38s
hybrid-pacbio-illumina | HG003  | 2h 15m
