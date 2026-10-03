CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Isrc
SRC := src/mps.cpp src/simplex.cpp src/presolve.cpp src/verify.cpp src/solve.cpp src/gpu.cpp
OBJ := $(SRC:.cpp=.o)

# Optional GPU offload of the pricing step:  make CUDA=1   (needs nvcc and a CUDA toolkit)
# --fmad=false keeps the GPU arithmetic identical to the CPU loop.
ifeq ($(CUDA),1)
NVCC ?= nvcc
OBJ := $(filter-out src/gpu.o,$(OBJ)) src/gpu_cpu_part.o src/gpu_cuda.o
CXXFLAGS += -DBS_WITH_CUDA
LDLIBS += -lcudart
endif

bin/bharatsolve: src/main.cpp $(OBJ)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp $(OBJ) $(LDLIBS)

%.o: %.cpp src/bs.hpp src/gpu.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifeq ($(CUDA),1)
src/gpu_cpu_part.o: src/gpu.cpp src/bs.hpp src/gpu.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
src/gpu_cuda.o: src/gpu_cuda.cu src/bs.hpp src/gpu.hpp
	$(NVCC) -O2 -std=c++17 --fmad=false -Isrc -c $< -o $@
endif

test: bin/bharatsolve tests/test_core
	./tests/test_core

test-slow: bin/bharatsolve tests/test_core
	BS_SLOW=1 ./tests/test_core

tests/test_core: tests/test_core.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ tests/test_core.cpp $(OBJ) $(LDLIBS)

# Tests the host side of the CUDA backend against a CPU stand-in for the CUDA runtime (no GPU needed).
CORE_OBJ := src/mps.o src/simplex.o src/presolve.o src/verify.o src/solve.o
test-cuda-mock: $(CORE_OBJ) tests/test_cuda_mock.cpp src/gpu_cuda.cu src/gpu.cpp
	$(CXX) $(CXXFLAGS) -DBS_WITH_CUDA -c src/gpu.cpp -o /tmp/bs_gpu_cpu_part.o
	$(CXX) $(CXXFLAGS) -x c++ -DBS_MOCK_CUDA -Itests/mock_cuda -c src/gpu_cuda.cu -o /tmp/bs_gpu_cuda_mock.o
	$(CXX) $(CXXFLAGS) -o /tmp/test_cuda_mock tests/test_cuda_mock.cpp $(CORE_OBJ) /tmp/bs_gpu_cpu_part.o /tmp/bs_gpu_cuda_mock.o
	/tmp/test_cuda_mock

clean:
	rm -f $(OBJ) src/gpu.o src/gpu_cpu_part.o src/gpu_cuda.o bin/bharatsolve tests/test_core
