CC := gcc
CXX := g++
FLAGS := -std=c++2a -O3 -march=native -Wall -Wextra -Wpedantic -fconcepts -fopenmp
LDFLAGS := -lyaml-cpp

GCCVERSION := $(shell expr `gcc -dumpversion | cut -f1 -d.`)

ifeq "${GCCVERSION}" "8"
	LDFLAGS += -lstdc++fs
endif

NW := dpu_sets
NW16S := dpu_16S

.PHONY: all clean 16s

all: ${NW} ${NW16S}

16s: ${NW16S}

clean:
	$(RM) ${NW}
	$(RM) ${NW16S}
	cd ./libnwdpu/dpu && make clean

SRC := ./src/main_sets.cpp ./src/fasta.cpp ./libnwdpu/host/dpu_common.cpp
SRC16S := ./src/main_16s.cpp ./src/fasta.cpp ./libnwdpu/host/dpu_common.cpp

${NW}: ${SRC}
	${CXX} ${FLAGS} $^ -o $@ ${LDFLAGS} `dpu-pkg-config --cflags --libs dpu`
	cd ./libnwdpu/dpu && make affine

${NW16S}: ${SRC16S}
	${CXX} ${FLAGS} $^ -o $@ ${LDFLAGS} `dpu-pkg-config --cflags --libs dpu`
	cd ./libnwdpu/dpu && make 16s
