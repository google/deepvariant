# DeepVariant training data

### WGS models

version | Replicates                             | #examples
------- | -------------------------------------- | -----------
v0.4    | 9 HG001                                | 85,323,867
v0.5    | 9 HG001<br>2 HG005<br>78 HG001 WES<br>1 HG005 WES<sup>[(1)](#vfootnote1)</sup> | 115,975,740
v0.6    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+     | 156,571,227
v0.7    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+     | 158,571,078
v0.8    | 12 HG001 PCR-free<br>2 HG005 PCR-free<br>4 HG001 PCR+<br>(and, more `dowsample_fraction` since last version)     | 346,505,686
v0.9    | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>2 HG006 PCR-free<br>2 HG007 PCR-free<br>5 HG001 PCR+     | 325,202,093
v0.10   | 10 HG001 PCR-free<br>2 HG005 PCR-free<br>2 HG006 PCR-free<br>2 HG007 PCR-free<br>5 HG001 PCR+     | 339,410,078
v1.0    | 11 HG001<br>2 HG005-HG007<br>2 HG002-HG004<sup>[(7)](#vfootnote7)</sup>     | 317,486,837
v1.1    | 12 HG001<br>3 HG002<br>3 HG004<br>3 HG005<br>3 HG006<br>3 HG007<sup>[(9)](#vfootnote9)</sup> | 388,337,190
v1.2    | 12 HG001<br>6 HG002<sup>[(12)](#vfootnote12)</sup><br>6 HG004<sup>[(12)](#vfootnote12)</sup><br>3 HG005<br>3 HG006<br>3 HG007 | 518,709,296
v1.3    | Same model as v1.2
v1.4    | 12 HG001<br>6 HG002<sup>[(12)](#vfootnote12)</sup><br>6 HG004<sup>[(12)](#vfootnote12)</sup><br>3 HG005<br>3 HG006<br>3 HG007 | 517,209,566
v1.5    | 13 HG001<br>14 HG002<br>8 HG004<br>9 HG005<br>4 HG006<br>4 HG007 | 815,200,320
v1.6    | 21 HG001<br>17 HG002<br>8 HG004<br>9 HG005<br>4 HG006<br>4 HG007 | 929,199,066

### WES models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v0.5    | 78 HG001<br>1 HG005 | 15,714,062
v0.6    | 78 HG001<br>1 HG005<sup>[(2)](#vfootnote2)</sup> | 15,705,449
v0.7    | 78 HG001<br>1 HG005 | 15,704,197
v0.8    | 78 HG001<br>1 HG005<sup>[(3)](#vfootnote3)</sup> | 18,683,247
v0.9    | 81 HG001<br>1 HG005<sup>[(3)](#vfootnote3)[(4)](#vfootnote4)[(5)](#vfootnote5)</sup> | 61,953,965
v0.10   | Same model as v0.9
v1.0    | 32 HG001<br>9 HG002<br>6 HG003<br>6 HG004<br>12 HG005<br>9 HG006<br>9 HG007<sup>[(7)](#vfootnote7)</sup> | 10,716,281
v1.1    | 41 HG001<br>9 HG002<br>6 HG004<br>12 HG005<br>9 HG006<br>9 HG007<sup>[(9)](#vfootnote9)</sup> | 13,450,688
v1.2    | 41 HG001<br>9 HG002<br>9 HG004<br>12 HG005<br>9 HG006<br>9 HG007<sup>[(11)](#vfootnote11)</sup> | 22,288,064
v1.3    | Same model as v1.2
v1.4    | 41 HG001<br>9 HG002<br>9 HG004<br>12 HG005<br>9 HG006<br>9 HG007<sup>[(11)](#vfootnote11)</sup> | 21,212,424
v1.5    | 40 HG001<br>9 HG002<br>9 HG004<br>12 HG005<br>9 HG006<br>9 HG007 | 21,027,625
v1.6    | 57 HG001<br>9 HG002<br>9 HG004<br>12 HG005<br>9 HG006<br>9 HG007 | 21,027,614

### PACBIO models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v0.8    | 16 HG002 | 160,025,931
v0.9    | 49 HG002 <sup>[(6)](#vfootnote6)</sup> | 357,507,235
v0.10   | 49 HG002, 2 HG003, 2 HG004, 1 HG002 (amplified) <sup>[(6)](#vfootnote6)</sup> | 472,711,858
v1.0    | 1 HG001<br>2 HG002<br>2 HG003<br>2 HG004<br>1 HG005 <sup>[(8)](#vfootnote8)</sup>  | 302,331,948
v1.1    | 1 HG001<br>9 HG002<br>2 HG004<br>1 HG005<sup>[(9)](#vfootnote9)</sup> | 569,225,616
v1.2    | 1 HG001<br>19 HG002<br>2 HG004<br>1 HG005<sup>[(10)](#vfootnote10)</sup> | 1,036,056,726
v1.3    | 1 HG001<br>19 HG002<br>3 HG004<br>1 HG005<br>1 HG006<br>1 HG007 | 1,177,109,190
v1.4    | 1 HG001<br>19 HG002<br>3 HG004<br>1 HG005<br>1 HG006<br>1 HG007 | 1,177,596,708
v1.5    | 3 HG001<br>29 HG002<br>7 HG004<br>2 HG005<br>3 HG006<br>2 HG007 | 1,729,659,396
v1.6    | 6 HG001<br>60 HG002<br>16 HG004<br>4 HG005<br>6 HG006<br>4 HG007 | 3,195,507,862

### ONT models

version | Replicates                  | #examples
------- | --------------------------- | ------------------------------
v1.6    | 3 HG001<br>1 HG004<br>1 HG005 | 534,302,654

### HYBRID models

version | Replicates                                               | #examples
------- | -------------------------------------------------------- | -----------
v1.0    | 10 HG002<br> 1 HG004<br> 1 HG005<br> 1 HG006<br> 1 HG007 | 193,076,623
v1.1    | Same model as v1.0                                       |
v1.2    | 10 HG002<br> 1 HG004<br> 1 HG005<br> 1 HG006<br> 1 HG007 | 214,302,681
v1.3    | Same model as v1.2                                       |
v1.4    | 10 HG002<br> 1 HG004<br> 1 HG005<br> 1 HG006<br> 1 HG007 | 215,863,645
v1.5    | 10 HG002<br> 1 HG004<br> 1 HG005<br> 1 HG006<br> 1 HG007 | 215,863,664
v1.6    | 10 HG002<br> 1 HG004<br> 1 HG005<br> 1 HG006<br> 1 HG007 | 215,353,081

<a name="vfootnote1">(1)</a>: In v0.5, we experimented with adding whole exome
sequencing data into training data. In v0.6, we took it out because it didn't
improve the WGS accuracy.

<a name="vfootnote2">(2)</a>: The training data are from the same replicates as
v0.5. The number of examples changed because of the update in
[haplotype_labeler](https://github.com/google/deepvariant/tree/r0.6/deepvariant/labeler/haplotype_labeler.py).

<a name="vfootnote3">(3)</a>: In v0.8, we used the
[Platinum Genomes Truthset](https://github.com/Illumina/PlatinumGenomes) to
create more training examples outside the GIAB confident regions.

<a name="vfootnote4">(4)</a>: Previously, we split train/tune by leaving 3 WES
for tuning. Starting from this release, we leave out chr1 and chr20 from
training, and use chr1 for tuning.

<a name="vfootnote5">(5)</a>: Starting from this version, we padded (100bps on
both sides) of the capture BED and used that for generating training examples.
We also added more `downsample_fraction`.

<a name="vfootnote6">(6)</a>: (Before v1.0) PacBio is the only one we currently
uses HG002 in training and tuning.

<a name="vfootnote7">(7)</a>: In v1.0, we train on HG002-HG004 for WGS as well,
but only using examples from the region of NIST truth confident region v4.2
subtracting v3.3.2.

<a name="vfootnote8">(8)</a>: In v1.0, PacBio training data contains training
examples with haplotag sorted images and unsorted images.

<a name="vfootnote9">(9)</a>: In v1.1, we exclude HG003 from training. And we
use all NIST truth confident regions for HG001-HG007 (except for HG003) for
training. We've always excluded chr20-22 from training.

<a name="vfootnote10">(10)</a>: In v1.2, we include new PacBio training data
from Sequel II, Chemistry 2.2.

<a name="vfootnote11">(11)</a>: Between v1.1 and v1.2, we fixed an issue where
make_examples can generate fewer class 0 (REF) training examples than before.
This is the reason for more training examples in v1.2 when number of samples
didn't increase.

<a name="vfootnote12">(12)</a>: In v1.2, we created BAM files with 100bp reads
and 125bp reads by trimming to augment the training data.

## Training data:

See "[An Extensive Sequence Dataset of Gold-Standard Samples for Benchmarking and Development](https://doi.org/10.1101/2020.12.11.422022)"
for a publicly available set of data we released. Data download information can
be found in the supplementary material.
