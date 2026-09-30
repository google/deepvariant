---
layout: post
title: "Training species-specific models with the genome of a single individual: Demonstration on Canola"
date:   2026-09-30
description: "In this blog we detail how to train a species-specific DeepVariant model using only a single individual when sequenced with both PacBio and Illumina technology, using Canola as a demonstration."
img: "assets/images/2026-09-30-canola/figure6-errors.png"
authors: ["satwood", "spaulraj", "jcao", "srennybyfield", "ssriram", "awcarroll", "pichuan"]
authors_by_org:
  - org: "Corteva Agriscience"
    members: ["satwood", "spaulraj", "jcao", "srennybyfield", "ssriram"]
  - org: "Google"
    members: ["awcarroll", "pichuan"]
---

## Summary

In this blog we detail how to train a species-specific DeepVariant model using only a single individual when sequenced with both PacBio and Illumina technology, using Canola as a demonstration. We quantify a ~40% error reduction for the re-trained model on the truth set. We show better concordance with orthogonal gene marker sets, and demonstrate improved variant discovery on an agriculturally interesting gene. The model is available for research use with DeepVariant.

## Introduction

The ability to sequence a species genome is essential to understand its traits. The methods to analyze the most standard sequence methods (mapping and variant detection to a reference genome) have been developed with mostly human genomes in mind, where gold standard datasets and resources (both data, existing tools, and funding incentive) is greatest. However, species differ substantially in their genome and population structure, including features like repeat or transposable element content, ploidy, divergence from the reference, among other factors.

Creating a species-optimized version of [DeepVariant](https://github.com/google/deepvariant) would enable strong performance regardless of genome characteristics. Previously, we demonstrated how a father-mother-child (trio) pair could be used to [generate a training set in mosquitoes](https://google.github.io/deepvariant/posts/2018-12-05-improved-non-human-variant-calling-using-species-specific-deepvariant-models/). External groups greatly expanded on this demonstration, including University of Otago training a [Kākāpō-specific model](https://www.nature.com/articles/s41559-023-02165-y), and investigators at University of Missouri [developing TrioTrain to train models for cattle](https://genome.cshlp.org/content/35/8/1859). Though these approaches are strong, they require a sequencing of a direct trio, which is not always available in projects.

Here, we demonstrate an effective strategy to train a DeepVariant model for a species without use of a trio by combining Pacific Biosciences and Illumina data for the same individual. We demonstrate this on Canola, an agriculturally important species for cooking oil production and protein-rich feed for animals.

## Methods - Overview

Here our goal is to train an Illumina DeepVariant model for Canola. Our input data will be the Illumina BAM files, but the key is to determine the Truth labels for training. For this, we will rely heavily on PacBio variant calls on the same individual.

Training a DeepVariant model requires creating a Truth VCF and a set of Confident BED regions. Our approach to this is to run the human release DeepVariant models (with the `--disable_small_model` flag) on both Illumina and PacBio BAM files. We take the PacBio VCF as the Truth variants, and then go through a process of excluding regions where we are not very confident they are correct.

To construct the confident regions, first, we make an inclusion BED file consisting of regions with 10 or more reads with MAPQ > 40 and 2 or fewer reads with MAPQ < 10.

Then we create an exclusion BED by writing a BED entry with +/- 50bp around any position that meets the following exclusion criteria:

* Any insertion or deletion of more than 50 bp (reference field in the VCF greater or equal to 50bp in length, or any Alternate allele greater or equal to 50bp in length) is excluded.
* For PacBio, any position with a no-call variant (genotype `./.`), that is a reference call with GQ < 20, is excluded.
* For PacBio, any alternate variant called with a GQ of < 30 is excluded.
* For SNPs, a variant position must have at least some evidence in the Illumina file, though the Illumina variant does not need to be called. So at least enough evidence to generate a candidate for DeepVariant in Illumina (that is more than 2 reads supporting a variant with an allele fraction > 0.12). If the variant does not have this level of evidence in Illumina, it is excluded.
* For Indels or multi-allelic variants with one or more alleles being Indel, both PacBio and Illumina must make the exact same variant call (position, reference, alternate, and genotype) call with PacBio calling at GQ > 30.

After creating the exclusion BED, the final Confident BED file is created by subtracting the exclusion BED file from the coverage inclusion BED file generated with the MAPQ statistics.

With these training VCF from the PacBio calls and the confident BED regions, training a DeepVariant model proceeds by following the established training workflow described by the general DeepVariant model. We train this model by warmstarting from the Human WGS DeepVariant model. All work is done with DeepVariant v1.9.

For all trainings, two chromosomes N17 and N18 were set aside from training data as the evaluation data used to select the best model. The chromosome N19 was fully set aside from all components of training as a hold out in every sample. The remaining chromosomes were used for training.

For this work, we used 5 individuals from different lines for training and fully held out 2 individuals from different lines as an evaluation set. Because the lines are proprietary to Corteva, we describe them with identifiers.

## Looking at Training Inputs and Preparation

To map Illumina samples, we used Bowtie2. PacBio sample mapping used pbmm2. The coverage of PacBio samples acquired was 20x–30x, which allows for a high confidence for the truth variants. The coverage of Illumina sequencing was 8x–15x, which reflects the coverages used for at-scale sequencing.

Mapping with Bowtie2 uses the following command:

```bash
bowtie2-align-s --wrapper basic-0 -p ${THREADS} -x ${REF} -X 2000 --no-discordant -1 ${READS} -2 ${MATES}
```

Mapping with pbmm2 uses the following command:

```bash
pbmm2 align ${REF} ${READS} ${BAM} --preset CCS --sort
```

![Figure 1]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure1-coverage.png)

*Figure 1. Coverage statistics for lines used in training and testing.*

As Figure 2 for zygosity indicates, these are inbred lines, and so overwhelmingly homozygous, though heterozygous sites do exist. We observe that even this small number of heterozygous sites still allow DeepVariant to call heterozygotes.

![Figure 2]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure2-zygosity.png)

*Figure 2. Number of variant calls in the PacBio sample by zygosity.*

The ultimate confident regions determined were consistent per sample, around 650 Mb, which is about 70% of the full size of the Canola genome.

![Figure 3]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure3-confident-regions-size.png)

*Figure 3. Size of the confident regions per line.*

![Figure 4]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure4-confident-variants.png)

*Figure 4. Number of confident variants per line.*

## Results - Improved accuracy on the test datasets

To assess the accuracy, we first used the ability to improve as measured against the truth labels derived from PacBio variant calls. Shown below are the precision and recall of the human DeepVariant v1.9 on Canola relative to the re-trained model. For SNPs, our Recall improves from ~0.94 to ~0.96 and our precision from ~0.92 to ~0.96. For Indels, our Recall improves from ~0.97 to ~0.98 and our precision improves from ~0.90 to ~0.98.

![Figure 5]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure5-precision-recall.png)

*Figure 5. Precision and Recall of the Illumina sequencing for the 2 held-out lines on the PacBio truth set (genome wide).*

We can also look at the total errors made on these lines and their origin:

![Figure 6]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure6-errors.png)

*Figure 6. The number of errors made per line on the Illumina sequencing of the 2 held-out lines relative to the PacBio Truth set. Errors are stratified by False Negative (FN) and False Positive (FP) and whether the error is a SNP or Indel.*

## Measuring Improved Accuracy on Corteva’s Pipeline and Markers

Direct comparisons between the re-trained DeepVariant model and an alternative variant pipeline that utilizes samtools pileup, referred here as SNPtool, demonstrated substantially improved recall while maintaining strong concordance with existing calls. The re-trained model recovered approximately 98% of all variant calls reported by SNPtool, which was consistent with expectation based on the highly conservative tuning of SNPtool for high recall. In addition, the re-trained model identified approximately 8% additional variants that were not detected by SNPtool.

To further compare variant calls between two methods with a fully orthogonal set, we used a well-curated set of Canola marker sites (referred here as GenoDB). Recall on calls observed at GenoDB marker sites improved from approximately 86% with SNPtool to approximately 92% with the re-trained DeepVariant model, suggesting that the variants newly discovered by the re-trained DeepVariant model are real variants which would be important to call.

![Figure 7]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure7-genodb.png)

*Figure 7 summarizes the overlap and divergence between variant calls generated by the two pipelines at GenoDB known marker positions, highlighting both the high concordance and improved recall from additional variants uniquely identified by the re-trained model.*

## Improved detection of heterozygous variants

Across evaluations, the re-trained DeepVariant model showed significantly improved performance in identifying heterozygous SNP and Indel calls, a known challenge for the internal SNPtool pipeline (due to specific tuning for homozygous calls) and for Canola due to genome polyploidy and variable sequencing depth. The PacBio-informed training strategy reduced misclassification of heterozygous sites and improved recall in regions that were previously difficult to resolve using short-read based DeepVariant re-training strategies. This improvement was observed consistently across benchmark datasets and internal validation analyses.

## Application to a Corteva Canola output trait use case

As an example of the practical importance of improving variants discovery, the re-trained DeepVariant model was applied to glucosinolate content (GLUC), an output trait use case of interest to Corteva, involving approximately 170 re-sequenced canola lines spanning coverage levels from ~3× to ~35×. Across this population-scale analysis, the re-trained model consistently identified more candidate variants per sample than SNPtool, particularly within a genomic loci of relevance for glucosinolate content. At the same time, the re-trained model preserved strong compatibility with existing marker sets that have previously been defined by pan-genomic analyses and marker validations, enabling seamless integration with downstream breeding workflows. As a specific illustration, the re-trained DV model identified SNPs at 64 out of 96 polymorphic sites defined by pan-genomic analysis, whereas SNPtool was only able to identify SNPs at 18 sites. Further, four sites that were previously vetted and validated through lab-based marker assays were all supported by SNP calls from re-trained DV, including at one of the four sites that was missed by SNPtool.

![Figure 8]({{ site.baseurl }}/assets/images/2026-09-30-canola/figure8-gluc.png)

*Figure 8. Comparison of variant call counts across GLUC samples generated using the re-trained DeepVariant model and SNPtool, illustrating the increased recall from discovery of expected homozygous and heterozygous candidate variants by the re-trained model enabled by the collaboration.*

## Downloading and Using the Model

We release the canola-retrained model used in this work for research use. Model files are available at the following bucket with public access and free download:

```
gs://brain-genomics-public/DeepVariant-blog/canola
```

This model was trained with DeepVariant v1.9, and should be used with that version of DeepVariant.

Please see this detailed walkthrough for how to download the model and how to run it using an example Canola dataset:

[https://gist.github.com/pichuan/c7e2ad4395343962a67fb7fa491a9d1f](https://gist.github.com/pichuan/c7e2ad4395343962a67fb7fa491a9d1f)

## Technical Notes

All analyses were executed at scale using Google Cloud compute resources, with observed runtimes consistent with expectations for population-scale variant calling (processing takes around 90 minutes on a 64-core machine). The re-trained model performed robustly across a wide range of sequencing depths; however, samples at the lowest coverage would benefit from post-calling filtering to balance sensitivity and specificity. In addition, historical PacBio datasets generated using older sequencing chemistry were expected to show reduced consistency relative to newer HiFi data, motivating consideration of re-sequencing with updated technology.

## Conclusion

We hope that the training recipe detailed here is straightforward to follow, and will add to the tools available to train species-specific models. As genomics enters a new era of scalability, our ability to understand the natural world to improve agriculture, medicine, and conservation will have new possibilities. Investigations such as this point to the fact that by investigating how to apply findings and methods across species, we can realize more of the potential insights in genomics data, and achieve more of the promise in this new scientific revolution.
