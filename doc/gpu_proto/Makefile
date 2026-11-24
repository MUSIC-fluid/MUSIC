# (2+1)D ideal KT prototype: CPU (OpenMP) and Metal backends on macOS,
# CUDA backend with nvcc (make cuda) on a machine that has it.
CXX      ?= clang++
OMP_DIR  ?= /opt/homebrew/opt/libomp
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wno-unused-parameter
OMPFLAGS  = -Xpreprocessor -fopenmp -I$(OMP_DIR)/include
OMPLIBS   = -L$(OMP_DIR)/lib -lomp
HERE     := $(abspath .)

all: kt2d

kt2d: kt2d_main.cpp kt2d_metal.mm kt2d.h kt2d_metal.h
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -DKT2D_WITH_METAL -c kt2d_main.cpp -o kt2d_main.o
	$(CXX) $(CXXFLAGS) -fobjc-arc -DKT2D_HEADER_PATH='"$(HERE)/kt2d.h"' -c kt2d_metal.mm -o kt2d_metal.o
	$(CXX) kt2d_main.o kt2d_metal.o -o $@ $(OMPLIBS) -framework Metal -framework Foundation

# type check of the CUDA source without an NVIDIA toolchain: clang's CUDA
# front end, host and device (sm_80) sides, with the runtime API stubbed
CUDA_CHECK = $(CXX) -x cuda -nocudainc -nocudalib --cuda-gpu-arch=sm_80 \
             -std=c++17 -Wall -DKT2D_NO_CUDA_SDK -fsyntax-only
cuda-syntax: kt2d_cuda.cu kt2d.h
	$(CUDA_CHECK) --cuda-host-only kt2d_cuda.cu
	$(CUDA_CHECK) --cuda-device-only kt2d_cuda.cu

cuda: kt2d_cuda.cu kt2d.h
	nvcc -O3 -std=c++17 -arch=native kt2d_cuda.cu -o kt2d_cuda

clean:
	rm -f *.o kt2d kt2d_cuda
