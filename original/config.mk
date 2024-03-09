# config include file for fletcher
# common for all backends

# Compilers
GCC=gcc
OMPSS2=/home/pyssouza/ompss2/build_llvm/bin/clang
#NVCC=nvcc
#PGCC=pgcc
#CLANG=clang

# Library paths
GCC_LIBS=-lm
#NVCC_LIBS=-lcudart -lstdc++    # it may include CUDA lib64 path...
#PGCC_LIBS=-lm
CLANG_LIBS=-lm

OMP_FLAG = -fopenmp

# PAPI flags
#PAPI_LIBS=-lpapi
