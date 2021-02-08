---
layout: post
title: "Improving Variant Calling using Haplotype Information"
date: 2021-02-08
description: "We discuss a new channel in DeepVariant which encodes haplotype information in long-read data, and was released with DeepVariant v1.1. We review how haplotypes relate to variant calling, show examples improved by the channel, and quantify the accuracy improvement with PacBio HiFi reads."
img: "assets/images/hp-channel/examples.png"
authors: ["danielecook","koles","pichuan","awcarroll"]
---

In this blog, we discuss a new channel in DeepVariant which encodes haplotype information in long-read data, and was released with DeepVariant v1.1. We review how haplotypes relate to variant calling, show examples improved by the channel, and quantify the accuracy improvement with PacBio HiFi reads.

## Introduction

Each individual inherits 23 pairs of chromosomes, one set from their mother and the other from their father. Genetic variation on each chromosome represents a haplotype, or group of “linked” alleles that are inherited from a single parent. Whether two parts of sequence are on the same haplotype (called _cis_) or opposite (_trans_) is important in [gene regulation](https://en.wikipedia.org/wiki/Cis-regulatory_element) and interpreting the [impact of variants](https://en.wikipedia.org/wiki/Compound_heterozygosity). 

In the process of sequencing, the genome is sampled in small pieces, which results in a loss of  long-range haplotype information. However, long-read sequencing approaches cover large enough regions to allow local reconstruction of which variants in the genome are on the same haplotype (in phase). Because human experts often use information about parental inheritance to interpret variants, we reasoned that this information would improve accuracy of DeepVariant.

DeepVariant classifies the genotype of candidate variant positions by representing attributes of the sequence data as a tensor.  The tensor is essentially a multidimensional image, consisting of channels that correspond to sequence features such as read base, base quality, and mapping quality. Each channel is centered on the candidate and has a width of 221 bp. This genomic interval provides enough local context to achieve high-accuracy with variant calling. 

<img src="{{ site.baseurl }}/assets/images/hp-channel/channels.png" alt="Examples of DeepVariant Channels" style="border: 1px solid black;">
<figcaption style='text-align: center;'>A subset of channels used by DeepVariant. Each channel is 221 bp wide.</figcaption>

The width of this window was originally determined for Illumina short reads. However, long-reads contain distant variant sites outside of the window that can be used to infer haplotypes (through a process called phasing). In the context of variant calling, haplotype information can provide evidence for or against a putative variant by linking similar evidence together, as opposed to random errors. Consider the example below of heterozygous variants flanking a homopolymorphic site where haplotype information provides support for a variant.

![Pileup Example CCS]({{ site.baseurl }}/assets/images/hp-channel/examples.png)
<figcaption>Alleles for variant 1 and 3 are segregated by parental haplotype. Alleles for variant 2 are ambiguous, but when examined in the context of haplotype we observe deletion alleles more often segregate with parent B reads. This additional information can be helpful to DeepVariant when classifying candidate variants.</figcaption>

Previously, we reported that sorting reads by haplotype reduced false negatives and positives by ~30% in _[Wenger 2019](https://www.nature.com/articles/s41587-019-0217-9)_ (and discussed in [this blog post](https://google.github.io/deepvariant/posts/2019-01-14-highly-accurate-snp-and-indel-calling-on-pacbio-ccs-with-deepvariant/)). These improvements were added to the [DeepVariant 1.0 release](https://ai.googleblog.com/2020/09/improving-accuracy-of-genomic-analysis.html). Notably, these improvements used an implicit representation of haplotype. We decided to investigate whether a haplotype channel, which provides a more direct representation of haplotype information, may further improve accuracy.

## Approach

### Integrating the Haplotype Channel

![DeepVariant CCS Workflow]({{ site.baseurl }}/assets/images/hp-channel/dv_ccs_calling.png)
<figcaption>Following alignment, DeepVariant calls heterozygous sites. Heterozygous sites are used by WhatsHap to perform phasing and haplotagging of reads. A second round of variant calling by DeepVariant incorporates the phase information from haplotagged reads.</figcaption>

To add a haplotype channel to DeepVariant, we perform two rounds of variant calling. The first round of variant calling identifies heterozygous variants. These heterozygous variants are used by [WhatsHap](https://whatshap.readthedocs.io/en/latest/guide.html) to perform phasing. WhatsHap will then annotate individual reads with an HP tag (a “haplotag”). The HP tag takes on values of 0, 1, or 2. `HP=1` and `HP=2` distinguish haplotypes at a particular locus, whereas `HP=0` indicates that no haplotype could be inferred. DeepVariant fetches the HP tag, and then constructs the haplotype channel.


<div style='text-align: center;'>
    <strong>chr1:2600683 12 bp deletion</strong>
</div>


<img src="{{ site.baseurl }}/assets/images/hp-channel/haplotype_channel.png" alt="The HP Channel" style="border: 1px solid black;">
<figcaption>A subset of channels. The haplotype channel is shown in the middle panel, with three groups of reads. The top (black) shows reads where no haplotype could be assigned (HP=0). The middle gray (HP=1) and bottom white (HP=2) reads represent different haplotypes. The deletion is only observed in the HP=2 haplotype.</figcaption>

### Training Models with the Haplotype Channel


In order to test whether a haplotype channel improves accuracy using PacBio HiFi reads, we generated the following training datasets from HG001, HG002, HG004, and HG005 Genome in a Bottle (GiaB) BAMs haplotagged using WhatsHap:

1.  **Baseline**: the existing model (no haplotype information included)
1.  **HP Sorted**: Reads sorted by haplotype (similar to _Wenger et al. 2019_)
1.  **HP Channel**: Reads sorted by haplotype and a haplotype channel.
1.  **Baseline & HP Channel**: Our baseline training set was combined with the HP Channel training set.

The first 3 experiments are designed to give a clear test of individual hypotheses. The last experiment is designed to compare with our production training process, as our release model needs to perform well on both haplotagged and untagged data. Having a single model that operates on both types of data also means only one model is needed to perform both rounds of variant calling when integrating haplotype information.

![HP-Channel Experiments]({{ site.baseurl }}/assets/images/hp-channel/experiments.png)
<figcaption>A visual overview of the experiments. Note that in the <i>Baseline & HP Channel</i> experiment, we add a blank channel to maintain compatibility between datasets. Compatibility between datasets allows us to combine training data with and without haplotype information, and develop a model that generalizes to both types of data.</figcaption>

## Results

![HP-Channel Experiments]({{ site.baseurl }}/assets/images/hp-channel/figure_1.png)
<figcaption style='text-align: center;'><strong>Figure 1</strong> Higher is better; F1 scores by SNP and INDEL across experimental conditions.</figcaption>

Overall, we observe the greatest improvements in INDEL calling (**Figure 1**). Sorting reads by haplotype drives significant improvement, but the addition of the haplotype channel further improves F1 scores.

![HP-Channel Experiments]({{ site.baseurl }}/assets/images/hp-channel/figure_2.png)
<figcaption style='text-align: center;'><strong>Figure 2</strong>  Lower is better; Relative improvement of DeepVariant PacBio models compared to the baseline model.</figcaption>

Another way to look at these results is to examine their relative improvement compared to the baseline. **Figure 2** shows the total number of errors (false positives + false negatives) identified when we compare DeepVariant with the truth set. Each bar lists the percent decrease compared to our baseline model. Here we can see more clearly that haplotype information does help reduce SNP calling errors. 

We reason that homopolymer INDEL errors are the main challenge in accurately calling variants in HiFi data. A variant caller has to distinguish the measurement variation in length (e.g. how many 7T or 8T runs) from the probability that one parent has 7Ts and the other 8Ts). By separating the signal by haplotype, the distributions become clearer. Labelling their boundaries with the explicit channel, and giving the variant caller clear information that the variant could be phased improves this further.

### Duplicate Training Comparison

One concern we had during the addition of the HP Channel was whether results would differ considerably between training runs.  To test this, we performed training and evaluation twice. The resulting evaluations did not differ more than 1.8% in terms of their total number of errors. 

### Runtime

We evaluated runtime using the [DeepVariant PacBio case study](https://github.com/google/deepvariant/blob/r1.0/docs/deepvariant-pacbio-model-case-study.md). The addition of the haplotype channel does not have a significant effect on runtime.

| Stage                            | Baseline (minutes)   | --use_hp_information  |
|:---------------------------------|:---------------------|:-----------------------|
| make_examples                    | 113m                 | 116m                   |
| call_variants                    | 218m                 | 222m                   |
| postprocess_variants (with gVCF) | 73m                  | 69m                    |
| **total**                            | **404m = 6.7 hours**     | **407m = 6.7 hours**       |

## Conclusion

Here we have demonstrated that the addition of a haplotype channel improves the accuracy of variant calling on PacBio HiFi data. These improvements have been integrated as part of the [DeepVariant 1.1 release](https://github.com/google/deepvariant/releases/tag/v1.1.0). See the [PacBio HiFi case study](https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-pacbio-model-case-study.md) for details on how to call variants using haplotype information.

This work is an example of how we can improve accuracy of analysis methods by using domain intuition about a problem and representing that information in a manner that a neural network can learn. Looking forward, we will investigate additional signals including mapping percent and base quality of long reads to see whether additional new channels lead to further improvements in accuracy. We look forward to continuing to improve the accuracy of DeepVariant for our users in future releases.
