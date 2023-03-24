CC := dpu-upmem-dpurte-clang
FLAGSP := -O3 -fno-builtin -DNR_TASKLETS=24 -DNR_GROUPS=6 -fshort-enums -DSTACK_SIZE_DEFAULT=384
FLAGS16S := -O3 -fno-builtin -DNR_TASKLETS=24 -DNR_GROUPS=6 -fshort-enums -DSTACK_SIZE_DEFAULT=384

NWP := nw_affine
NW16S := nw_16s

.PHONY: all affine clean 16s

all:${NWP} ${NW16S}

affine: ${NWP}

16s: ${NW16S}

clean:
	$(RM) ${NWP}
	$(RM) ${NW16S}

SRCP := ${NWP}.c
SRC16S := ${NW16S}.c
	
${NWP}: ${SRCP}
	${CC} ${FLAGSP} $^ -o $@

${NW16S}: ${SRC16S}
	${CC} ${FLAGS16S} $^ -o $@