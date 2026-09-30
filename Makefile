CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Isrc
SRC := src/mps.cpp src/simplex.cpp src/presolve.cpp src/verify.cpp src/solve.cpp
OBJ := $(SRC:.cpp=.o)

bin/bharatsolve: src/main.cpp $(OBJ)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp $(OBJ)

%.o: %.cpp src/bs.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: bin/bharatsolve tests/test_core
	./tests/test_core

tests/test_core: tests/test_core.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ tests/test_core.cpp $(OBJ)

clean:
	rm -f $(OBJ) bin/bharatsolve tests/test_core
