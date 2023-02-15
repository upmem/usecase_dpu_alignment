# Needleman and Wunsch on DPU

## Libraries needed

> - UPMEM SDK
> - c++17 with filesystem library (gcc vesrion >=8)
> - yaml-cpp

## Build

> make

## Run application

> ./dpu_alignment

Alignment parameters can be changed in `params.yaml`
Output scores and cigars in scores.txt and cigars.txt respectively.


## Dataset format

must have the form:

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