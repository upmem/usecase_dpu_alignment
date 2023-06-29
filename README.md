# Needleman and Wunsch on DPU

This repository implements two different alignment workload on UPMEM PiM.

The first one named dpu_alignment compares sets of sequences with all sequences in the same set being compared with all the other in the same set.

The second one named dpu_16s is an exemple of comparing a database of 16S RNA, all against all comparison. This can be use for more than just 16S RNA sequences.

## Libraries needed

> - UPMEM SDK
> - filesystem and span std library (gcc vesrion >=10)
> - yaml-cpp

## Build

> make

## Run application

For set comparison:
> ./dpu_alignment

For 16S RNA:
> ./dpu_16s

Alignment parameters can be changed in `params.yaml` and `16s.yaml`
Output scores and cigars in scores.txt and cigars.txt respectively.


## Dataset format

### Set comparison fasta file form

```
> set {set_number} ...
{sequence}
> set {set_number} ...
{sequence}
.
.
.
```
All sequences with same number will be pair-aligned, same set number sequences must be contiguous !

### 16S comparison fasta file form:

```
> ...
{sequence}
> ...
{sequence}
.
.
.
```
Limited to 1024 sequences for now !
