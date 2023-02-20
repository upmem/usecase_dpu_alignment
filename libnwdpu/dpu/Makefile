CC := dpu-upmem-dpurte-clang
FLAGSP := -O3 -fno-builtin -DNR_TASKLETS=24 -DNR_GROUPS=6 -fshort-enums -DSTACK_SIZE_DEFAULT=384

NWP := nw_affine

.PHONY: all affine clean

all:${NWP}

affine: ${NWP}

clean:
	$(RM) ${NWP}

SRCP := ${NWP}.c
	
${NWP}: ${SRCP}
	${CC} ${FLAGSP} $^ -o $@