# Build

## libraries

> - UPMEM SDK
> - c++17 filesystem library
> - yaml-cpp

Makefile

> make

# Getting Started

run  `./dpu_alignment`

Alignment parameters can be changed in `params.yaml`


# dataset

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

All sequence with same number will be pair-aligned, same set number sequences must be contiguous !