---
layout: post
title: "Improving DeepVariant Accuracy by Providing Chemistry Priors"
date:   2026-10-26
description: "In this blog we demonstrate an improvement to DeepVariant accuracy across PacBio chemistries by combining chemistry metadata with a Feature-wise Linear Modulation (FiLM) addition to the DeepVariant model."
img: "assets/images/2026-10-26-chemistry-priors/figure2-f1.png"
authors: ["kmuret", "danielecook", "pichuan", "awcarroll", "koles"]
---

## Introduction

Deep learning models learn from the underlying distribution of their input data how to accurately perform the tasks for which they are trained. [DeepVariant](https://google.github.io/deepvariant/posts/2020-02-20-looking-through-deepvariants-eyes/) uses a Convolutional Neural Network to perform the task of identifying whether, and how, a sample’s genetic sequence differs from a reference genome. Doing this accurately requires understanding the error profile of the sequencing instrument: knowing what noise to expect is what lets the model distinguish a true variant from an artifact.

Training a single variant calling model across diverse datasets requires balancing broad generalization with adaptation to version-specific error patterns. PacBio HiFi sequencing has progressed through several generations, including [Sequel II (2019)](https://www.pacb.com/press_releases/pacific-biosciences-launches-new-sequel-ii-system-featuring-8-times-the-dna-sequencing-data-output/), [Revio v1 (2023)](https://www.pacb.com/press_releases/pacbio-announces-revio-a-revolutionary-new-long-read-sequencing-system-designed-to-provide-15-times-more-hifi-data-and-human-genomes-at-scale-for-under-1000/), and the newer [SPRQ (2024)](https://www.pacb.com/press_releases/pacbio-announces-sprq-chemistry-for-revio-sequencing-systems-a-major-advance-reducing-the-cost-of-a-hifi-human-genome-to-less-than-500/) and [SPRQ-Nx (2026)](https://www.pacb.com/blog/improving-hifi-sequencing-accuracy-with-google-deepconsensus-and-alphaevolve/) chemistries. Although all of these platforms produce high-accuracy reads with similar characteristics, each generation changes aspects of the polymerase kinetics, sequencing kits, and DNA input requirements. These changes shift the empirical error distributions, particularly in homopolymers (consecutive runs of the same base), which are the most challenging context for PacBio sequencing and where errors typically appear as insertions and deletions (INDELs).

We simultaneously want DeepVariant to perform well across all chemistries without burdening the user with choosing among chemistry-specific models. Historically, handling multiple chemistries required either training separate models for each chemistry or combining all datasets into a single unconditioned model, which is the approach that has been taken for the DeepVariant PacBio models to date. Combining all datasets into an unconditioned model allows the network to learn shared features from a larger pool of examples, but forces a single set of weights to average across different error modes. In this work, we investigate whether conditioning a single DeepVariant model on sequencing chemistry metadata using Feature-wise Linear Modulation (FiLM) can capture chemistry-specific error profiles while sharing a unified model across all PacBio datasets.

## Methods

### Integrating Chemistry Metadata with FiLM

In the standard DeepVariant architecture, the InceptionV3 backbone outputs a 2048-dimensional feature vector (<i>h</i><sub>img</sub>) after global average pooling, which passes through a dropout layer and a linear classification head. To modulate the extracted image features directly, we implemented Feature-wise Linear Modulation (FiLM) between the global average pooling layer and the classification head (Figure 1). The model takes two categorical metadata inputs: the sequencing chemistry and a binary PCR-amplification flag. Each is mapped through its own 32-dimensional embedding layer, and the two embeddings are concatenated and batch-normalized to produce a 64-dimensional metadata vector <i>m</i>. In our experiments the amplification flag produced no measurable improvement (evaluation results were unchanged with the flag on or off), so this work focuses on chemistry. Two linear layers project <i>m</i> into 2048-dimensional scale (<i>γ</i>(<i>m</i>)) and shift (<i>β</i>(<i>m</i>)) vectors, which modulate the pooled image representation element-wise before classification:

<p style="text-align: center;"><i>h</i><sub>FiLM</sub> = <i>γ</i>(<i>m</i>) ⊙ <i>h</i><sub>img</sub> + <i>β</i>(<i>m</i>)</p>

![Figure 1]({{ site.baseurl }}/assets/images/2026-10-26-chemistry-priors/figure1-architecture.png)

*Figure 1. Comparison of the standard unconditioned DeepVariant architecture (left) and the chemistry-conditioned FiLM architecture (right).*

### Data

During candidate generation in `make_examples`, DeepVariant parses the `BindingKit`, `SequencingKit`, and `CM` fields from the `@RG` lines of the BAM header and maps each dataset to a chemistry category (Table 1).

| Group | Binding Kit | Sequencing Kit | CM Tag / String |
| --- | --- | --- | --- |
| `SEQUEL_II` | 102-194-100 | 101-826-100 | S/P5-C2/... |
| `REVIO_V1` | 102-739-100 | 102-118-800 | R/P1-C1/... |
| `REVIO_SPRQ` | 103-496-900 | 103-496-700 | R/P2-C2/... |
| `REVIO_SPRQ_Nx` | 103-496-900 | 103-726-100 | R/P2-C3/... |

*Table 1. Mapping of PacBio BAM header fields to the chemistry categories used in this work. Kit part numbers shown are examples; each chemistry has several.*

To test whether chemistry conditioning improves accuracy, we trained two models on the same set of 24 Genome in a Bottle (GIAB) PacBio HiFi BAMs, with chromosomes 20–22 held out from training: an unconditioned Baseline model and a chemistry-conditioned FiLM model (Exp FiLM). The HG002 Sequel II, Revio v1, and SPRQ datasets used for evaluation below were part of this training set. For SPRQ-Nx, we used three public PacBio HG002 datasets (one SMRT Cell used three times, released in 2026) that were not used in training. We also evaluated on HG003 (Sequel II, Revio v1, and SPRQ), a sample completely held out from training.

## Results

We evaluated both models on whole-genome HG002 data across four PacBio chemistries: Sequel II, Revio v1, SPRQ, and SPRQ-Nx. To compare chemistries without confounding differences in depth, each dataset was downsampled to ~30x coverage (Revio v1 was used at its full ~28x) and evaluated against the v5 GIAB Truth Set (Figure 2).

![Figure 2]({{ site.baseurl }}/assets/images/2026-10-26-chemistry-priors/figure2-f1.png)

*Figure 2. Whole-genome INDEL (left) and SNP (right) F1 scores on HG002 at ~30x coverage across four PacBio chemistries (three SPRQ-Nx datasets, Revio v1 at its full ~28x) for the unconditioned Baseline model and the Exp FiLM model. Higher is better.*

Conditioning on chemistry metadata improves INDEL F1 scores over the Baseline model across all four chemistries, with the largest INDEL gains on SPRQ, followed by Sequel II and SPRQ-Nx, while SNP F1 scores are essentially unchanged (+0.01 percentage points). To understand whether these F1 improvements come from fewer false positives or fewer false negatives, we examined precision and recall for each chemistry at ~30x coverage (Figure 3).

![Figure 3]({{ site.baseurl }}/assets/images/2026-10-26-chemistry-priors/figure3-precision-recall.png)

*Figure 3. Precision and recall comparison between the Baseline model (open circles) and the Exp FiLM model (filled circles) at ~30x coverage (Revio v1 at its full ~28x). SPRQ-Nx is shown for one of the three datasets (Use1).*

Across all four chemistries, the Exp FiLM model shifts performance toward higher precision, with recall higher or essentially unchanged. Rather than trading one metric for the other, conditioning on chemistry reduces false positives without meaningfully increasing false negatives.

Looking at the total number of variant calling errors (false positives plus false negatives) across the whole genome shows the scale of these reductions (Figure 4). At ~30x coverage, the Exp FiLM model reduces total errors across all four chemistries compared to the Baseline model. Since homopolymer errors are the main source of difference between generations, providing the chemistry label yields the largest error reductions on INDELs (3.1% to 31.2% fewer errors), with smaller reductions in SNP errors on HG002 (2.7% to 6.0% fewer).

![Figure 4]({{ site.baseurl }}/assets/images/2026-10-26-chemistry-priors/figure4-errors.png)

*Figure 4. Total whole-genome INDEL errors (left) and SNP errors (right) at ~30x coverage (Revio v1 at its full ~28x) comparing the Baseline model and the Exp FiLM model. Lower is better. SPRQ-Nx is shown for one of the three datasets (Use1).*

To see how chemistry conditioning performs at lower sequencing depths, we evaluated both models across downsampled coverages ranging from 35x down to 9x (Revio v1: 28x to 7x) (Table 2). On SPRQ, Sequel II, and SPRQ-Nx, the Exp FiLM model maintains higher INDEL F1 scores and fewer total errors than the Baseline model across the entire coverage range, with gains for SPRQ and Sequel II generally growing as coverage decreases, while SPRQ-Nx gains stay small (+0.12 to +0.26 percentage points). On Revio v1, the Exp FiLM model reduces total errors down to 15x, while the unconditioned Baseline model has a slightly higher INDEL F1 at 15x and below.

<style>
.t2 { border: 1px solid #dadce0; border-radius: 8px; overflow: hidden; margin: 8px 0 12px; font-size: 14px; }
.t2-head { display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 4px 16px; background: #f8f9fa; padding: 12px 14px; }
.t2-title { font-weight: 700; color: #202124; font-size: 15px; }
.t2-sub { color: #5f6368; font-size: 12px; }
.t2-scroll { overflow-x: auto; }
.t2 table { margin: 0; width: 100%; border: 0; border-collapse: collapse; color: #202124; text-align: left; }
.t2 table tr:nth-child(n) { background: #fff; }
.t2 table tr.t2-alt { background: #fafafa; }
.t2 table th, .t2 table td { border: 0; padding: 5px 8px; white-space: nowrap; vertical-align: middle; }
.t2 table th { background: #f1f3f4; color: #3c4043; font-size: 11px; font-weight: 700; letter-spacing: 0.04em; text-transform: uppercase; line-height: 1.3; }
.t2 table tr.t2-first td { border-top: 1px solid #bdc1c6; }
.t2 table td.t2-chem { background: #fff; font-weight: 700; vertical-align: top; padding-top: 7px; }
.t2-dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 6px; }
.t2-cov, .t2-base, .t2-exp, .t2-pill { font-family: "Roboto Mono", Menlo, Consolas, monospace; font-size: 13px; }
.t2-cov, .t2-exp { font-weight: 700; }
.t2-base { color: #5f6368; }
.t2-pill { display: inline-block; min-width: 64px; padding: 1px 6px; border-radius: 4px; text-align: center; font-weight: 700; }
.t2-good { background: #e6f4ea; color: #137333; }
.t2-bad { background: #fce8e6; color: #c5221f; }
.t2-saved { background: #e8f0fe; color: #1a73e8; }
.t2-note { border-top: 1px solid #bdc1c6; color: #5f6368; font-size: 12px; line-height: 1.5; padding: 10px 14px; }
</style>
<div class="t2">
<div class="t2-head"><span class="t2-title">Coverage Titration: Variant F1-Score &amp; Accuracy Performance Matrix</span><span class="t2-sub">GIAB HG002 WGS &bull; Baseline vs. Exp FiLM</span></div>
<div class="t2-scroll">
<table>
<thead><tr><th>Chemistry</th><th>Coverage<br>(depth)</th><th>Baseline<br>INDEL F1</th><th>Exp FiLM<br>INDEL F1</th><th>INDEL<br>&Delta;F1 (pp)</th><th>SNP<br>&Delta;F1 (pp)</th><th>Errors saved<br>(INDEL + SNP)</th></tr></thead>
<tbody>
<tr class="t2-first"><td class="t2-chem" rowspan="6" style="color: #ea4335;"><span class="t2-dot" style="background: #ea4335;"></span>Sequel II</td><td class="t2-cov">35x (1.00)</td><td class="t2-base">0.969771</td><td class="t2-exp">0.971699</td><td><span class="t2-pill t2-good">+0.19</span></td><td><span class="t2-pill t2-good">+0.007</span></td><td><span class="t2-pill t2-saved">+4,054</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">29x (0.85)</td><td class="t2-base">0.963292</td><td class="t2-exp">0.966334</td><td><span class="t2-pill t2-good">+0.30</span></td><td><span class="t2-pill t2-good">+0.009</span></td><td><span class="t2-pill t2-saved">+6,124</span></td></tr>
<tr><td class="t2-cov">24x (0.70)</td><td class="t2-base">0.952642</td><td class="t2-exp">0.957318</td><td><span class="t2-pill t2-good">+0.47</span></td><td><span class="t2-pill t2-good">+0.015</span></td><td><span class="t2-pill t2-saved">+9,514</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">19x (0.55)</td><td class="t2-base">0.934013</td><td class="t2-exp">0.941403</td><td><span class="t2-pill t2-good">+0.74</span></td><td><span class="t2-pill t2-good">+0.018</span></td><td><span class="t2-pill t2-saved">+14,510</span></td></tr>
<tr><td class="t2-cov">14x (0.40)</td><td class="t2-base">0.896660</td><td class="t2-exp">0.908168</td><td><span class="t2-pill t2-good">+1.15</span></td><td><span class="t2-pill t2-good">+0.026</span></td><td><span class="t2-pill t2-saved">+22,161</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">9x (0.25)</td><td class="t2-base">0.807032</td><td class="t2-exp">0.821229</td><td><span class="t2-pill t2-good">+1.42</span></td><td><span class="t2-pill t2-good">+0.025</span></td><td><span class="t2-pill t2-saved">+25,206</span></td></tr>
<tr class="t2-first"><td class="t2-chem" rowspan="6" style="color: #1a73e8;"><span class="t2-dot" style="background: #1a73e8;"></span>Revio v1</td><td class="t2-cov">28x (1.00)</td><td class="t2-base">0.965044</td><td class="t2-exp">0.966111</td><td><span class="t2-pill t2-good">+0.11</span></td><td><span class="t2-pill t2-good">+0.009</span></td><td><span class="t2-pill t2-saved">+2,619</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">24x (0.85)</td><td class="t2-base">0.957007</td><td class="t2-exp">0.957819</td><td><span class="t2-pill t2-good">+0.08</span></td><td><span class="t2-pill t2-good">+0.014</span></td><td><span class="t2-pill t2-saved">+2,522</span></td></tr>
<tr><td class="t2-cov">20x (0.70)</td><td class="t2-base">0.943746</td><td class="t2-exp">0.944357</td><td><span class="t2-pill t2-good">+0.06</span></td><td><span class="t2-pill t2-good">+0.015</span></td><td><span class="t2-pill t2-saved">+2,230</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">15x (0.55)</td><td class="t2-base">0.920894</td><td class="t2-exp">0.920760</td><td><span class="t2-pill t2-bad">-0.01</span></td><td><span class="t2-pill t2-good">+0.013</span></td><td><span class="t2-pill t2-saved">+809</span></td></tr>
<tr><td class="t2-cov">11x (0.40)</td><td class="t2-base">0.872803</td><td class="t2-exp">0.872104</td><td><span class="t2-pill t2-bad">-0.07</span></td><td><span class="t2-pill t2-good">+0.008</span></td><td><span class="t2-pill t2-bad">-712</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">7x (0.25)</td><td class="t2-base">0.754643</td><td class="t2-exp">0.753676</td><td><span class="t2-pill t2-bad">-0.10</span></td><td><span class="t2-pill t2-good">+0.006</span></td><td><span class="t2-pill t2-bad">-2,423</span></td></tr>
<tr class="t2-first"><td class="t2-chem" rowspan="6" style="color: #f29900;"><span class="t2-dot" style="background: #f29900;"></span>SPRQ</td><td class="t2-cov">37x (1.00)</td><td class="t2-base">0.977056</td><td class="t2-exp">0.984011</td><td><span class="t2-pill t2-good">+0.70</span></td><td><span class="t2-pill t2-good">+0.013</span></td><td><span class="t2-pill t2-saved">+13,337</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">31x (0.85)</td><td class="t2-base">0.970845</td><td class="t2-exp">0.979920</td><td><span class="t2-pill t2-good">+0.91</span></td><td><span class="t2-pill t2-good">+0.014</span></td><td><span class="t2-pill t2-saved">+17,102</span></td></tr>
<tr><td class="t2-cov">26x (0.70)</td><td class="t2-base">0.961958</td><td class="t2-exp">0.973182</td><td><span class="t2-pill t2-good">+1.12</span></td><td><span class="t2-pill t2-good">+0.016</span></td><td><span class="t2-pill t2-saved">+21,038</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">20x (0.55)</td><td class="t2-base">0.946970</td><td class="t2-exp">0.959864</td><td><span class="t2-pill t2-good">+1.29</span></td><td><span class="t2-pill t2-good">+0.018</span></td><td><span class="t2-pill t2-saved">+24,002</span></td></tr>
<tr><td class="t2-cov">15x (0.40)</td><td class="t2-base">0.916283</td><td class="t2-exp">0.931006</td><td><span class="t2-pill t2-good">+1.47</span></td><td><span class="t2-pill t2-good">+0.017</span></td><td><span class="t2-pill t2-saved">+26,699</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">9x (0.25)</td><td class="t2-base">0.838603</td><td class="t2-exp">0.851979</td><td><span class="t2-pill t2-good">+1.34</span></td><td><span class="t2-pill t2-good">+0.011</span></td><td><span class="t2-pill t2-saved">+21,904</span></td></tr>
<tr class="t2-first"><td class="t2-chem" rowspan="6" style="color: #0d652d;"><span class="t2-dot" style="background: #0d652d;"></span>SPRQ-Nx</td><td class="t2-cov">35x (0.76)</td><td class="t2-base">0.981848</td><td class="t2-exp">0.983660</td><td><span class="t2-pill t2-good">+0.18</span></td><td><span class="t2-pill t2-good">+0.012</span></td><td><span class="t2-pill t2-saved">+4,094</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">30x (0.65)</td><td class="t2-base">0.977522</td><td class="t2-exp">0.979556</td><td><span class="t2-pill t2-good">+0.20</span></td><td><span class="t2-pill t2-good">+0.013</span></td><td><span class="t2-pill t2-saved">+4,513</span></td></tr>
<tr><td class="t2-cov">25x (0.53)</td><td class="t2-base">0.970048</td><td class="t2-exp">0.972485</td><td><span class="t2-pill t2-good">+0.24</span></td><td><span class="t2-pill t2-good">+0.012</span></td><td><span class="t2-pill t2-saved">+5,105</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">19x (0.42)</td><td class="t2-base">0.956163</td><td class="t2-exp">0.958772</td><td><span class="t2-pill t2-good">+0.26</span></td><td><span class="t2-pill t2-good">+0.017</span></td><td><span class="t2-pill t2-saved">+5,715</span></td></tr>
<tr><td class="t2-cov">14x (0.30)</td><td class="t2-base">0.926550</td><td class="t2-exp">0.928932</td><td><span class="t2-pill t2-good">+0.24</span></td><td><span class="t2-pill t2-good">+0.015</span></td><td><span class="t2-pill t2-saved">+4,922</span></td></tr>
<tr class="t2-alt"><td class="t2-cov">9x (0.19)</td><td class="t2-base">0.841332</td><td class="t2-exp">0.842523</td><td><span class="t2-pill t2-good">+0.12</span></td><td><span class="t2-pill t2-good">+0.006</span></td><td><span class="t2-pill t2-saved">+922</span></td></tr>
</tbody>
</table>
</div>
<div class="t2-note">Note: Benchmarked on whole-genome GIAB HG002 v5.0q benchmark set. &Delta;F1 is in percentage points (pp). Positive &Delta;F1 and positive Errors Saved denote improvements over Baseline.<br>SPRQ-Nx: public PacBio 2026Q2 HG002 data (not used in training), downsampled so the top row is ~35x.</div>
</div>

*Table 2. Comparison of variant calling accuracy and total errors saved between the Baseline model and the Exp FiLM model across downsampled coverages from 35x to 9x on HG002 (Revio v1: 28x to 7x). SPRQ-Nx results are shown for one of the three datasets (Use1).*

Because chromosomes 20–22 were held out of training, we also evaluated HG002 on those chromosomes alone at ~30x. Summed across Sequel II, Revio v1, and SPRQ, the Exp FiLM model made 9.7% fewer INDEL errors on this held-out region (11,087 to 10,012), with reductions on every chemistry (Sequel II -4.8%, Revio v1 -0.8%, SPRQ -26.0%). SNP errors were slightly higher (+3.5%; 7,017 to 7,263), with small increases on every chemistry (Sequel II +2.9%, Revio v1 +3.6%, SPRQ +4.3%). For comparison, across the whole genome the same three runs had 13.1% fewer INDEL errors and 3.6% fewer SNP errors.

As a stronger test, we also evaluated both models on HG003, a sample completely held out from training, at full coverage (Sequel II 34x, Revio v1 32x) on the GIAB v4.2.1 benchmark (Sequel II, Revio v1, and SPRQ; no HG003 SPRQ-Nx data was available). Summed across the three runs, the Exp FiLM model made 13.6% fewer INDEL errors (20,883 to 18,037), with reductions on every chemistry (Sequel II -1.2%, Revio v1 -22.1%, SPRQ -13.4%). SNP errors were essentially unchanged overall (+0.1%; 20,059 to 20,085), higher on Revio v1 (+7.9%) and lower on Sequel II (-2.9%) and SPRQ (-2.2%). Because the HG003 benchmark (GIAB v4.2.1) covers fewer difficult regions than the HG002 v5 benchmark, the two sets of results are not directly comparable.

Finally, to confirm that these improvements come from the chemistry label itself rather than the extra parameters in the FiLM layer, we ran a negative control experiment on the HG002 Revio v1 dataset (Figure 5). At inference time, we compared the Baseline model against the Exp FiLM model given the correct Revio v1 label and the same Exp FiLM model given an incorrect SPRQ label. While the true chemistry label reduces both INDEL and SNP errors compared to the Baseline, giving the model the wrong chemistry label causes INDEL errors to increase substantially while SNP errors return to slightly above Baseline levels (+3.5%). This indicates that the model uses the chemistry label to adjust its INDEL calls to match the error profile of the specified chemistry.

![Figure 5]({{ site.baseurl }}/assets/images/2026-10-26-chemistry-priors/figure5-negative-control.png)

*Figure 5. Whole-genome INDEL errors (left) and SNP errors (right) on HG002 Revio v1 data at ~30x coverage (same data as Figure 4) for the Baseline model (gray), the Exp FiLM model with the correct Revio v1 label (blue), and the Exp FiLM model with an incorrect SPRQ label (red). Lower is better.*

## Conclusion

Conditioning DeepVariant on PacBio sequencing chemistry with Feature-wise Linear Modulation (FiLM) improves INDEL accuracy across Sequel II, Revio v1, SPRQ, and SPRQ-Nx datasets within a single unified model, including SPRQ-Nx data not used in training. In early experiments, adding the chemistry label only just before the final classification layer (late fusion) gave little improvement over the Baseline, which suggests that the label needs to change how the image features themselves are used, as FiLM does. By reading kit identifiers directly from BAM `@RG` headers during candidate generation, the model adapts to the error profile of each chemistry without splitting the training data across separate chemistry-specific checkpoints. This approach simplifies model maintenance as new sequencing kits are released and provides a general way to incorporate additional run-level metadata into DeepVariant, such as the instrument software version, which can also change error modes. We plan to include this chemistry-aware model in a future DeepVariant release.

## Acknowledgements

We thank Aaron Wenger from PacBio for information about PacBio run metadata, and for several helpful comments and suggestions in drafting this blog.
