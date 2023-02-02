CC := gcc
CXX := g++
FLAGS := -std=c++2a -O3 -march=native -Wall -Wextra -Wpedantic -fconcepts -fopenmp
LDFLAGS := -lyaml-cpp

GCCVERSION := $(shell expr `gcc -dumpversion | cut -f1 -d.`)

ifeq "${GCCVERSION}" "8"
	LDFLAGS += -lstdc++fs
endif

NW := dpu_alignment

.PHONY: all clean

all: ${NW}

clean:
	$(RM) ${NW}
	cd ./libnwdpu/dpu && make clean

SRC := ./src/main.cpp ./src/fasta.cpp ./libnwdpu/host/dpu_common.cpp

${NW}: ${SRC}
	${CXX} ${FLAGS} $^ -o $@ ${LDFLAGS} `dpu-pkg-config --cflags --libs dpu`
	cd ./libnwdpu/dpu && make clean && make affine