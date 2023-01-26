CC := gcc
CXX := g++
FLAGS := -std=c++2a -O3 -march=native -Wall -Wextra -Wpedantic -fconcepts -fopenmp

NW := dpu_alignment

.PHONY: all clean

all: ${NW}

clean:
	$(RM) ${NW}
	cd ./libnwdpu/dpu && make clean

SRC := ./src/main.cpp ./src/fasta.cpp ./libnwdpu/host/dpu_common.cpp

${NW}: ${SRC}
	${CXX} ${FLAGS} $^ -o $@ -lyaml-cpp `dpu-pkg-config --cflags --libs dpu`
	cd ./libnwdpu/dpu && make clean && make affine