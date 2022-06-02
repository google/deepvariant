# DeepTrio training data

### WGS models

| version      | Replicates                         | #examples   |
| ------------ | ---------------------------------- | ----------- |
| Child model  |                                    |             |
| 1.1.0        | 4 HG001/NA12891/NA12892 trios<br>7 HG005/HG006/HG007 trios <br>3 HG002/HG003/HG004 trios| 566,589,652<sup>[(1)](#vfootnote1)</sup> |
| 1.2.0        | (Same model as 1.1.0)              |             |
| 1.3.0        | (Same model as 1.1.0)              |             |
| 1.4.0        | 4 HG001/NA12891/NA12892 trios<br>7 HG005/HG006/HG007 trios <br>3 HG002/HG003/HG004 trios | 704,228,446 |
| Parent model |                                    |             |
| 1.1.0        | 7 HG005/HG006/HG007 trios <br> 3 HG002/HG003/HG004 trios | 315,847,934 |
| 1.2.0        | (Same model as 1.1.0)              |             |
| 1.3.0        | (Same model as 1.1.0)              |             |
| 1.4.0        | 7 HG005/HG006/HG007 trios <br> 3 HG002/HG003/HG004 trios | 457,374,516 |

### WES models

| version      | Replicates                                      | #examples  |
| ------------ | ----------------------------------------------- | ---------- |
| Child model  |                                                 |            |
| 1.1.0        | 27 HG001/NA12891/NA12892 trios<br>6 HG005/HG006/HG007 trios <br>7 HG002/HG003/HG004 trios  | 18,002,596 |
| 1.2.0        | (Same model as 1.1.0)              |             |
| 1.3.0        | (Same model as 1.1.0)              |             |
| 1.4.0        | 27 HG001/NA12891/NA12892 trios<br>6 HG005/HG006/HG007 trios <br>6 HG002/HG003/HG004 trios  | 27,776,416 |
| Parent model |                                                 |            |
| 1.1.0        | 6 HG005/HG006/HG007 trios <br> 6 HG002/HG003/HG004 trios  | 4,131,018  |
| 1.2.0        | (Same model as 1.1.0)              |             |
| 1.3.0        | (Same model as 1.1.0)              |             |
| 1.4.0        | 6 HG005/HG006/HG007 trios <br> 6 HG002/HG003/HG004 trios  | 13,036,995  |

### PACBIO models<sup>[(2)](#vfootnote2)</sup><sup>[(3)](#vfootnote3)</sup>

| version      | Replicates                         | #examples   |
| ------------ | ---------------------------------- | ----------- |
| Child model  |                                    |             |
| 1.1.0        | 1 HG005/HG006/HG007 trio <br>8 HG002/HG003/HG004 trios | 397,610,700 |
| 1.2.0        | 1 HG005/HG006/HG007 trio <br>8 HG002/HG003/HG004 trios | 406,893,180<sup>[(4)](#vfootnote4)</sup> |
| 1.3.0        | 2 HG005/HG006/HG007 trio <br>10 HG002/HG003/HG004 trios | 539,382,124<sup>[(5)](#vfootnote5)</sup> |
| 1.4.0        | (Same model as 1.3.0)              |             |
| Parent model |                                    |             |
| 1.1.0        | 1 HG005/HG006/HG007 trio <br> 8 HG002/HG003/HG004 trios | 386,418,918 |
| 1.2.0        | 1 HG005/HG006/HG007 trio <br>8 HG002/HG003/HG004 trios | 392,749,204<sup>[(4)](#vfootnote4)</sup> |
| 1.3.0        | 2 HG005/HG006/HG007 trio <br>10 HG002/HG003/HG004 trios | 533,353,050<sup>[(5)](#vfootnote5)</sup> |
| 1.4.0        | (Same model as 1.3.0)              |             |


<a name="vfootnote1">(1)</a>: We include HG002/HG003/HG004 for training WGS
model, but only using examples from the region of NIST truth confident region
v4.2 subtracting v3.3.2.

<a name="vfootnote2">(2)</a>: We use the entire HG002/HG003/HG004 trio for
PacBio model training.

<a name="vfootnote3">(3)</a>: PacBio training data contains training examples
with haplotag sorted images and unsorted images.

<a name="vfootnote4">(4)</a>: In v1.2.0, we updated the NIST truth versions we
used for training.

<a name="vfootnote5">(5)</a>: In v1.3.0, we included PacBio Sequel II Chemistry
v2.2 data in the training dataset. And we updated to NIST truth version to
v4.2.1.
