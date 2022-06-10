---
layout: post
title: "Adding Custom Channels to DeepVariant"
date: 2022-06-09
description: "This tutorial will show you how to add custom channels to DeepVariant, train new models, and use them to perform variant calling. Custom channels can be used to tailor DeepVariant to a specific sequencing platform or application."
img: "assets/images/custom-channels/channel_construction.png"
authors: ["lucasbrambrink", "danielecook"]
---

## Background

DeepVariant is a deep-learning based variant caller that formulates variant calling as an image classification problem. Aligned sequences are read from a [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file, and DeepVariant extracts a series of sequence features from read pileups and converts them to a corresponding set of input images that we call channels. Each channel encodes a sequence feature extracted from reads. These channels are the fundamental unit by which we feed data into our model, and we combine a series of them as input to train our model to accurately predict genotypes.

All DeepVariant production models use these [6 core channels](https://google.github.io/deepvariant/posts/2020-02-20-looking-through-deepvariants-eyes/) ([Figure 1](#fig1)):

* read base
* base quality
* mapping quality
* strand of alignment
* read supports variant
* base differs from ref

These get concatenated into a single image and are used as input into DeepVariant models.

<img id='fig1' src="{{ site.baseurl }}/assets/images/custom-channels/core_channels.png" alt="Examples of DeepVariant Channels" >
<figcaption style='text-align: center;'>Figure 1: DeepVariant's 6 core channels.</figcaption>

New sequencing technologies may produce novel signals (_e.g._ raw basecalling metrics) that could be encoded as a channel to improve accuracy specific to that platform. Alternatively, custom channels can specify derived properties (_e.g._ GC content) that can be experimented with to improve model accuracy on existing platforms.

To allow models to capture additional signals, we introduced a framework for developing custom channels in DeepVariant 1.2. Custom channels are specified as a list using the `--channels` flag with `make_examples` (see the original [Nature Biotechnology paper](https://www.nature.com/articles/nbt.4235) for an introduction to how DeepVariant works). Each custom channel can encode novel signals or derived sequence properties that might help DeepVariant make more accurate predictions. This blog post will cover how to develop a new custom channel, how to train a model for that channel, and how to run inference using a model that incorporates a custom channel.

We have successfully used the new framework with the release of DeepVariant v1.4 to include a new custom channel called `insert_size` which encodes the fragment length of Illumina paired end sequencing reads. Our analysis finds that inclusion of the `insert_size` channel in our models results in models that are 4.7% more accurate on WGS and 9.8% on WES.

In another recent development, Ultima Genomics has [put together custom channels](https://cdn.sanity.io/files/l7780ks7/production/385224b3861b445e96665a4ad212f8742746dd30.pdf) for their [new sequencing platform](https://www.biorxiv.org/content/10.1101/2022.05.29.493900v3). The Ultima Genomics sequencing platform produces probabilities for error in homopolymer length (encoded in the base-quality string, together with an auxiliary field called TP) and a probability score for missed signals (T0). These signals are incorporated as custom channels and lead to improved variant calling accuracy with their platform.

This tutorial will show you how to add custom channels to DeepVariant, train new models, and use them to perform variant calling. Users are expected to have intermediate experience with C++ and Python.

## How are Channels Constructed?

<img id='fig2' src="{{ site.baseurl }}/assets/images/custom-channels/channel_construction.png" alt="Examples of DeepVariant Channels" >
<figcaption style='text-align: center;'>Figure2: Channel Elements: Channels encode a reference strip, and read pileup.</figcaption>


Although we visualize channels as images, their underlying representation is as multidimensional tensor objects. Each channel is encoded as a matrix of values ranging from 0-255. Channels have two components ([Figure 2](#fig2)). The top five pixels are known as the “reference strip”, and can be used to specify information derived from the reference, although only a few channels use this component (_e.g._ read base). Below the reference strip is the read pileup. Given that channels are encoded using pileup information, the width corresponds to a set number of bases whereas the height reflects the maximum number of reads. While these parameters are configurable, experimentation has yielded the optimal width and height for short paired end sequence reads to be [221 and 100](https://github.com/google/deepvariant/blob/r1.3/deepvariant/dv_constants.py#L36-L40), respectively. `make_examples` generates these matrices, and stacks them together to form 3D tensor objects as model input.

To build a custom channel, we’ll have to define some functions and logic to extract the desired sequence features using the [Nucleus library](https://github.com/google/nucleus), and scale their values appropriately within a 0-255 value range. There are two channel types currently:

* __read-level channels__ - assign a single value to an entire read within the pileup (e.g. read supports variant).
* __base-level channels__ - base-level channels specify values at the base level (e.g. base quality).

## How to Add a New Custom Channel

The process of adding an Optchannel will require forking DeepVariant, allowing you to make edits to the codebase. The following example will add a new channel called `read_orientation` that encodes information about the orientation of paired end read alignments. Once the channel has been added, we will describe how to train a model and apply it to call variants using this new channel.

```bash
git clone https://github.com/google/deepvariant.git
cd deepvariant
```

__1. Declare a new channel option__

`deepvariant/dv_constants.py`

The `OPT_CHANNELS` list defines the set of available custom channels that can be used by make_examples.py, so we will begin by naming our channel and adding it to this list.

<div class="language-python highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="c1"># Define available OptChannels (optional extra channels).
</span><span class="n">OPT_CHANNELS</span> <span class="o">=</span> <span class="p">[</span>
    <span class="s">'read_mapping_percent'</span><span class="p">,</span> <span class="s">'avg_base_quality'</span><span class="p">,</span> <span class="s">'identity'</span><span class="p">,</span>
    <span class="s">'gap_compressed_identity'</span><span class="p">,</span> <span class="s">'gc_content'</span><span class="p">,</span> <span class="s">'is_homopolymer'</span><span class="p">,</span>
    <span class="s">'homopolymer_weighted'</span><span class="p">,</span> <span class="s">'blank'</span><span class="p">,</span> <span class="s">'insert_size'</span><span class="p">,</span><span class="highlight-code s">'read_orientation'</span><span class="p">,</span>
<span class="p">]</span>
</code></pre></div></div>

__2. Add this channel name to the DeepVariantChannelEnum__

`deepvariant/protos/deepvariant.proto`

This proto is used to help the python and C++ components of DeepVariant communicate with one another. By convention, enum values are prefixed with CH_ and uppercased. We’ll add an enum value to represent the new read orientation channel. Custom user channels should use enum values above 1000 so that it does not collide with future channel additions by the DeepVariant team.

<div class="language-proto highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="kd">enum</span> <span class="n">DeepVariantChannelEnum</span> <span class="p">{</span>
  <span class="c1">// Default should be unspecified.
</span>  <span class="na">CH_UNSPECIFIED</span> <span class="o">=</span> <span class="mi">0</span><span class="p">;</span>

  <span class="c1">// 6 channels that exist in all DeepVariant production models.
</span>  <span class="na">CH_READ_BASE</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span>
  <span class="na">CH_BASE_QUALITY</span> <span class="o">=</span> <span class="mi">2</span><span class="p">;</span>
  <span class="na">CH_MAPPING_QUALITY</span> <span class="o">=</span> <span class="mi">3</span><span class="p">;</span>
  <span class="na">CH_STRAND</span> <span class="o">=</span> <span class="mi">4</span><span class="p">;</span>
  <span class="na">CH_READ_SUPPORTS_VARIANT</span> <span class="o">=</span> <span class="mi">5</span><span class="p">;</span>
  <span class="na">CH_BASE_DIFFERS_FROM_REF</span> <span class="o">=</span> <span class="mi">6</span><span class="p">;</span>

  <span class="err">…</span>

  <span class="c1">// The following channels correspond to the "Opt Channels" defined in
</span>  <span class="c1">// deepvariant/pileup_channel_lib.h:
</span>  <span class="na">CH_READ_MAPPING_PERCENT</span> <span class="o">=</span> <span class="mi">11</span><span class="p">;</span>
  <span class="na">CH_AVG_BASE_QUALITY</span> <span class="o">=</span> <span class="mi">12</span><span class="p">;</span>
  <span class="na">CH_IDENTITY</span> <span class="o">=</span> <span class="mi">13</span><span class="p">;</span>
  <span class="na">CH_GAP_COMPRESSED_IDENTITY</span> <span class="o">=</span> <span class="mi">14</span><span class="p">;</span>
  <span class="na">CH_GC_CONTENT</span> <span class="o">=</span> <span class="mi">15</span><span class="p">;</span>
  <span class="na">CH_IS_HOMEOPOLYMER</span> <span class="o">=</span> <span class="mi">16</span><span class="p">;</span>
  <span class="na">CH_HOMEOPOLYMER_WEIGHTED</span> <span class="o">=</span> <span class="mi">17</span><span class="p">;</span>
  <span class="na">CH_BLANK</span> <span class="o">=</span> <span class="mi">18</span><span class="p">;</span>
  <span class="na">CH_INSERT_SIZE</span> <span class="o">=</span> <span class="mi">19</span><span class="p">;</span>
  <span class='highlight-code'><span class="na">CH_READ_ORIENTATION</span> <span class="o">=</span> <span class="mi">1000</span><span class="p">;</span></span>
<span class="p">}</span>
</code></pre></div></div>


__3. Add the channel name as a constant in our OptChannel library file__

`deepvariant/pileup_channel_lib.h`

This `pileup_channel_lib.h` file contains the logic for drawing custom channels. We will first define a value for the `read_orientation` channel to help bridge the gap between Python and C++.

<div class="language-cpp highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="c1">//--------------//</span>
<span class="c1">// Opt Channels //</span>
<span class="c1">//--------------//</span>

<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_read_mapping_percent</span> <span class="o">=</span> <span class="s">"read_mapping_percent"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_avg_base_quality</span> <span class="o">=</span> <span class="s">"avg_base_quality"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_identity</span> <span class="o">=</span> <span class="s">"identity"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_gap_compressed_identity</span> <span class="o">=</span> <span class="s">"gap_compressed_identity"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_gc_content</span> <span class="o">=</span> <span class="s">"gc_content"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_is_homopolymer</span> <span class="o">=</span> <span class="s">"is_homopolymer"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_homopolymer_weighted</span> <span class="o">=</span> <span class="s">"homopolymer_weighted"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_blank</span> <span class="o">=</span> <span class="s">"blank"</span><span class="p">;</span>
<span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_insert_size</span> <span class="o">=</span> <span class="s">"insert_size"</span><span class="p">;</span>
<span class='highlight-code'><span class="k">static</span> <span class="k">const</span> <span class="k">auto</span><span class="o">&amp;</span> <span class="n">ch_read_orientation</span> <span class="o">=</span> <span class="s">"read_orientation"</span><span class="p">;</span></span>
</code></pre></div></div>

__4. Map this channel name to the DeepVariantChannelProto enum__

`deepvariant/pileup_channel_lib.h`

A simple function is used to map the channel name to its corresponding enum. So we will follow the pattern of existing channels in this function for the `read_orientation` channel.

<div class="language-cpp highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="kr">inline</span> <span class="n">DeepVariantChannelEnum</span> <span class="nf">ChannelStrToEnum</span><span class="p">(</span><span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&amp;</span> <span class="n">channel</span><span class="p">)</span> <span class="p">{</span>
  <span class="c1">// Maps channel string representation to DeepVariantChannelEnum</span>
  <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_mapping_percent</span><span class="p">)</span>
    <span class="k">return</span> <span class="n">DeepVariantChannelEnum</span><span class="o">::</span><span class="n">CH_READ_MAPPING_PERCENT</span><span class="p">;</span>
  <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_avg_base_quality</span><span class="p">)</span>
    <span class="k">return</span> <span class="n">DeepVariantChannelEnum</span><span class="o">::</span><span class="n">CH_AVG_BASE_QUALITY</span><span class="p">;</span>
  <span class="p">...</span>
  <span class='highlight-code'><span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_orientation</span><span class="p">)</span>
    <span class="k">return</span> <span class="n">DeepVariantChannelEnum</span><span class="o">::</span><span class="n">CH_READ_ORIENTATION</span><span class="p">;</span></span>
  <span class="n">CHECK</span><span class="p">(</span><span class="nb">false</span><span class="p">)</span> <span class="o">&lt;&lt;</span> <span class="s">"Channel '"</span> <span class="o">&lt;&lt;</span> <span class="n">channel</span> <span class="o">&lt;&lt;</span> <span class="s">"' should have a corresponding "</span>
      <span class="o">&lt;&lt;</span> <span class="s">"enum in DeepVariantChannelEnum."</span><span class="p">;</span>
<span class="p">}</span>

</code></pre></div></div>

__5. Define the channel encoding__

`deepvariant/pileup_channel_lib.h`

This is naturally the most important step as it specifies how to encode your channel. Whatever data you are interested in representing, here is where you map it to a pixel value from 0 to 255. In our case, we will take 4 bitwise flags present in the `FLAG` section of a given read in a BAM file, and map them to a single pixel value per read.

```cpp
// Encodes read orientation and alignment information of read
inline std::vector<uint8> ReadOrientation(const Read& read) {
  /* Maps the bit fields `proper_placement`, `secondary_alignment`,
   * `supplementary_alignment` and `reverse_strand` into a 4-bit
   * integer value, i.e 0-15:
   *
   *   0000000
   *   0001000 <- proper_placement (8)
   *   0010000 <- secondary_alignment (16)
   *   0100000 <- supplementary_alignment (32)
   *   1000000 <- alignment.position.reverse_strand (64)
   *
   */
  int value = 0;
  if (read.proper_placement()) {
    value += 8;
  }
  if (read.secondary_alignment()) {
    value += 16;
  }
  if (read.supplementary_alignment()) {
    value += 32;
  }
  if (read.alignment().position().reverse_strand()) {
    value += 64;
  }
  std::vector<uint8> read_orientation_vector(read.aligned_sequence().size(),
                                             value);
  return read_orientation_vector;
}
```

__6. Tell OptChannels.CalculateChannels to use your new Channel__
`deepvariant/pileup_channel_lib.h`


<div class="language-cpp highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="k">class</span> <span class="nc">OptChannels</span> <span class="p">{</span>
 <span class="k">public</span><span class="o">:</span>
  <span class="n">std</span><span class="o">::</span><span class="n">map</span><span class="o">&lt;</span><span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="p">,</span> <span class="n">std</span><span class="o">::</span><span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">char</span><span class="o">&gt;&gt;</span> <span class="n">data_</span><span class="p">;</span>
  <span class="n">std</span><span class="o">::</span><span class="n">map</span><span class="o">&lt;</span><span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="p">,</span> <span class="n">std</span><span class="o">::</span><span class="n">vector</span><span class="o">&lt;</span><span class="kt">unsigned</span> <span class="kt">char</span><span class="o">&gt;&gt;</span> <span class="n">ref_data_</span><span class="p">;</span>
  <span class="kt">void</span> <span class="nf">CalculateChannels</span><span class="p">(</span><span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">vector</span><span class="o">&lt;</span><span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&gt;&amp;</span> <span class="n">channels</span><span class="p">,</span>
                         <span class="k">const</span> <span class="n">Read</span><span class="o">&amp;</span> <span class="n">read</span><span class="p">)</span> <span class="p">{</span>
    <span class="c1">// Calculates values for each channel</span>
    <span class="k">for</span> <span class="p">(</span><span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&amp;</span> <span class="n">channel</span> <span class="o">:</span> <span class="n">channels</span><span class="p">)</span> <span class="p">{</span>
      <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_mapping_percent</span><span class="p">)</span> <span class="p">{</span>
        <span class="n">data_</span><span class="p">[</span><span class="n">channel</span><span class="p">].</span><span class="n">assign</span><span class="p">(</span>
            <span class="p">{</span><span class="n">ScaleColor</span><span class="p">(</span><span class="n">ReadMappingPercent</span><span class="p">(</span><span class="n">read</span><span class="p">),</span> <span class="n">MaxMappingPercent</span><span class="p">)});</span>
      <span class="c1">// ...</span>
      <span class="p">}</span> <span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_insert_size</span><span class="p">)</span> <span class="p">{</span>
        <span class="n">data_</span><span class="p">[</span><span class="n">channel</span><span class="p">]</span> <span class="o">=</span> <span class="n">ReadInsertSize</span><span class="p">(</span><span class="n">read</span><span class="p">);</span>
      <span class="p">}</span>
      <span class='highlight-code'>
      <span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_orientation</span><span class="p">)</span> <span class="p">{</span>
        <span class="n">data_</span><span class="p">[</span><span class="n">channel</span><span class="p">]</span> <span class="o">=</span> <span class="n">ReadOrientation</span><span class="p">(</span><span class="n">read</span><span class="p">);</span></span>
      <span class="p">}</span>
    <span class="p">}</span>
  <span class="p">}</span>
<span class="err">…</span>
<span class="p">}</span>
</code></pre></div></div>

__7. Define how to represent the reference rows__

`deepvariant/pileup_channel_lib.h`

Each pileup image has 5 reference rows at the top. You can define how these are represented in your channel in CalculateRefRows. For most custom channels we simply use 5 blank (255) rows.

<div class="language-cpp highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="kt">void</span> <span class="nf">CalculateRefRows</span><span class="p">(</span><span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">vector</span><span class="o">&lt;</span><span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&gt;&amp;</span> <span class="n">channels</span><span class="p">,</span>
                      <span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&amp;</span> <span class="n">ref_bases</span><span class="p">)</span> <span class="p">{</span>
  <span class="c1">// Calculates reference row values for each channel</span>
  <span class="c1">// Create a fake read to represent reference bases.</span>
  <span class="n">Read</span> <span class="n">refRead</span><span class="p">;</span>
  <span class="k">for</span> <span class="p">(</span><span class="k">const</span> <span class="n">std</span><span class="o">::</span><span class="n">string</span><span class="o">&amp;</span> <span class="n">channel</span> <span class="o">:</span> <span class="n">channels</span><span class="p">)</span> <span class="p">{</span>
    <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_mapping_percent</span><span class="p">)</span> <span class="p">{</span>
      <span class="n">ref_data_</span><span class="p">[</span><span class="n">channel</span><span class="p">].</span><span class="n">assign</span><span class="p">({</span><span class="k">static_cast</span><span class="o">&lt;</span><span class="n">uint8</span><span class="o">&gt;</span><span class="p">(</span><span class="n">kMaxPixelValueAsFloat</span><span class="p">)});</span>
    <span class="c1">// ...</span>
    <span class='highlight-code'><span class="p">}</span> <span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_read_orientation</span><span class="p">)</span> <span class="p">{</span>
      <span class="n">ref_data_</span><span class="p">[</span><span class="n">channel</span><span class="p">].</span><span class="n">assign</span><span class="p">({</span><span class="k">static_cast</span><span class="o">&lt;</span><span class="n">uint8</span><span class="o">&gt;</span><span class="p">(</span><span class="n">kMaxPixelValueAsFloat</span><span class="p">)});</span></span>
    <span class="p">}</span> <span class="k">else</span> <span class="k">if</span> <span class="p">(</span><span class="n">channel</span> <span class="o">==</span> <span class="n">ch_gc_content</span><span class="p">)</span> <span class="p">{</span>
    <span class="c1">// ...</span>
    <span class="p">}</span>
  <span class="p">}</span>
<span class="p">}</span>
</code></pre></div></div>

That’s it! Now let’s see it in action.

## Verifying a Channel was Added Correctly
Now that we have added code for generating our new channel, we can build a docker image to test it out.

```bash
DOCKER_IMAGE="deepvariant/custom-channel"
docker build -t ${DOCKER_IMAGE} .
```

It is a good idea to start by checking that DeepVariant can use our new channel during make_examples. Assuming we have prepared the same environment as the WES case study, we can declare the following variables:

```bash
CHANNEL="read_orientation"
REF=/reference/GRCh38_no_alt_analysis_set.fasta
BAM=/input/HG003.novaseq.wes_idt.100x.dedup.bam
BED=/input/idt_capture_novogene.grch38.bed
EXAMPLES=/output/HG002.chr20.${CHANNEL}.examples.tfrecord.gz
```

Then we may run `make_examples`.

```bash
docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  ${DOCKER_IMAGE} \
  /opt/deepvariant/bin/make_examples \
    --mode calling \
    --ref ${REF} \
    --reads ${BAM} \
    --regions "chr20:0-500000" \
    --examples ${EXAMPLES} \
    --channels ${CHANNEL} \
    --logtostderr
```

Assuming that DeepVariant found some candidates and created examples for them, we can use show_examples to visualize them and ensure our new channel is part of the pileup image.

```bash
docker run \
  -v "${PWD}/output":"/output" \
  ${DOCKER_IMAGE_TAG}  \
  /opt/deepvariant/bin/show_examples \
    --examples="${EXAMPLES}" \
    --output="/output/" \
    --num_records=10 \
    --scale=10 \
    --image_type channels
```

The output of `show_examples` rasterizes each encoded variant example, and clearly shows the new we added (Figure 3). It worked!

<img id='fig3' src="{{ site.baseurl }}/assets/images/custom-channels/new_channel.png" alt="A New DeepVariant Channel" >
<figcaption style='text-align: center;'>Figure 3: A new Channel is visualized.</figcaption>

## Training with a Custom Channel

Now that we have determined that our new channel works, we can use it to train a new model. We will use data from the [WES case study](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-exome-case-study.md#prepare-environment), but closely follow the detailed training instructions of the [advanced training case study](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-training-case-study.md), highlighting here the most important aspects that pertain to the inclusion of our channel.

```bash
TRUTH_VCF="/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
```

First, we run `make_examples` again in `training` mode, which will properly label the examples based on the truth set and confident regions we provide. This allows us to generate a [training](https://developers.google.com/machine-learning/glossary/#training-set) and [validation set](https://developers.google.com/machine-learning/glossary/#validation-set), using different chromosomal regions.

```bash
docker run \
  -v "${PWD}/input":"/input" \
  -v "${PWD}/output":"/output" \
  -v "${PWD}/reference":"/reference" \
  -v "${PWD}/benchmark":"/benchmark" \
  ${DOCKER_IMAGE}  \
  /opt/deepvariant/bin/make_examples \
  --mode training \
  --ref ${REF} \
  --reads ${BAM} \
  --truth_variants "${TRUTH_VCF}" \
  --confident_regions "${TRUTH_BED}" \
  --examples "output//training_set.with_label.tfrecord@.gz" \
  --channels ${CHANNEL} \
  --regions "'chr1'" \
  --logtostderr
  ```

Following the same scheme of [shuffling our examples](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-training-case-study.md#shuffle-each-set-of-examples-and-generate-a-data-configuration-file-for-each), we produce a `dataset_config.pbtxt` file for each training and validation set.

```bash
TRAINING_DIR=/output/training_dir
DATASET_CONFIG_PBTXT="/output/training_set.dataset_config.pbtxt"
mkdir -p ${PWD}${TRAINING_DIR}
```

Now `model_train` can use this training data with our new custom channel and train its deep neural net.

<div class="language-bash highlighter-rouge"><div class="highlight"><pre class="highlight"><code>docker run <span class="se">\</span>
  <span class="nt">-v</span> <span class="s2">"</span><span class="k">${</span><span class="nv">PWD</span><span class="k">}</span><span class="s2">/output"</span>:<span class="s2">"/output"</span> <span class="se">\</span>
  <span class="k">${</span><span class="nv">DOCKER_IMAGE</span><span class="k">}</span> <span class="se">\</span>
  /opt/deepvariant/bin/model_train <span class="se">\</span>
  <span class="nt">--dataset_config_pbtxt</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">DATASET_CONFIG_PBTXT</span><span class="k">}</span><span class="s2">"</span> <span class="se">\</span>
  <span class="nt">--train_dir</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">TRAINING_DIR</span><span class="k">}</span><span class="s2">"</span> <span class="se">\</span>
  <span class="nt">--model_name</span><span class="o">=</span><span class="s2">"</span><span class="k">${</span><span class="nv">CHANNEL</span><span class="k">}</span><span class="s2">_model"</span> <span class="se">\</span>
  <span class="nt">--number_of_steps</span><span class="o">=</span>50000 <span class="se">\</span>
  <span class="nt">--save_interval_secs</span><span class="o">=</span>300 <span class="se">\</span>
  <span class="nt">--batch_size</span><span class="o">=</span>32 <span class="se">\</span>
  <span class='highlight-code'><span class="nt">--learning_rate</span><span class="o">=</span>0.0005
</span></code></pre></div></div>

It is possible to set the flags [--start_from_checkpoint](https://github.com/google/deepvariant/blob/r1.4/deepvariant/model_train.py#L106-L111) and [--allow_warmstart_from_different_num_channels](https://github.com/google/deepvariant/blob/r1.4/deepvariant/modeling.py#L107-L110) to start from an existing model checkpoint trained on a different channel configuration. This may significantly reduce the amount of additional training required.

As outlined in our [training guide](https://github.com/google/deepvariant/blob/r1.4/docs/deepvariant-training-case-study.md#start-model_train-and-model_eval), we can run `model_eval` simultaneously to determine the best model checkpoint. Once the training is done, we make sure to select our best performing model checkpoint.

```bash
MODEL_CHECKPOINT_PATH=$(cat "${TRAINING_DIR}"/best_checkpoint.txt)
```

## Running a Model with a New Custom Channel

With our newly trained model ready, we can now run DeepVariant to perform variant calling. Since we excluded chromosome 20 in our training set and validation sets, we can use it here to test the model’s accuracy.

```bash
docker run \
  -v "${PWD}/output":"/output" \
  ${DOCKER_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type WES \
  --customized_model "${MODEL_CHECKPOINT_PATH}" \
  --make_examples_extra_args=channels="${CHANNEL}" \
  --ref "${REF}" \
  --reads "${BAM}" \
  --regions "chr20" \
  --output_vcf "/output/test_set.vcf.gz"
```

Finally we may run the accuracy analysis as outlined in the [WES case study](https://github.com/google/deepvariant/blob/r1.3/docs/deepvariant-exome-case-study.md#benchmark-on-all-chromosomes) using [hap.py](https://github.com/Illumina/hap.py). If the data encoded in our new channel provides pertinent sequence information, we may reasonably expect some improvement in DeepVariant’s ability to correctly call variants.

## Conclusion

Custom channels are a powerful way to further improve DeepVariant’s variant calling performance and tailor models for specific sequencing platforms or applications. We hope you will be able to deploy custom channels to further improve DeepVariant performance.
