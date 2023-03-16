---
layout: post
title: "Variational inference for polygenic risk scores"
date:   2023-03-16
description: "We introduce a new way to build polygenic risk scores from GWAS summary
statistics using black box variational inference."
img: "assets/images/PacBio-2019-01/figure2b.png"
authors: ["nickfurlotte"]
---

On the Genomics team in Google Health AI, we think a lot about [polygenic risk scores (PRS)](https://www.genome.gov/Health/Genomics-and-Medicine/Polygenic-risk-scores). One concept we have been interested in is [black box variational inference (bbVI)](https://arxiv.org/abs/1401.0118) for crafting better PRS. Black box variational inference is a form of variational inference (VI) that solves the optimization problem using stochastic optimization and automatic differentiation. This is a cool technique because it allows the practitioner a high degree of flexibility on the statistical modeling side as compared to methods such as Markov-chain Monte Carlo (MCMC). If you agree, or even if you are only vaguely familiar with some words in this paragraph but think they sound enticing, please keep reading!

Given our interest in VI, we were excited to see a [preprint](https://www.biorxiv.org/content/10.1101/2022.05.10.491396v2) from Zabad et al. posted in May of 2022 that introduces the key concepts for using VI to create PRS. They provide great baseline results using a traditional VI approach. However, traditional VI does not employ the black box technique which translates into less modeling flexibility. We believe that bbVI implemented using libraries such as [TensorFlow Probability](https://www.tensorflow.org/probability) has a lot of potential to be a powerful and flexible methodology for building Bayesian methods, so we figured a blog post laying out some of the basics might be useful for the community. 

Here, we first give background context on PRS, LDpred, MCMC and variational inference. If you are already familiar with these concepts, please feel free to skip those sections and head right to the “bbVI + PRS = bbviPRS” section, where we introduce bbviPRS: a method to construct a PRS from GWAS summary statistics that uses black box variational inference.

