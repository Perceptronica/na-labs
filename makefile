CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra
INCLUDES = -I.

EXECUTABLES = lu_dec triagonal iterations jacobi

all: $(EXECUTABLES)

lu: lu-dec.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

triagonal: triagonal.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

iterations: iterations.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

jacobi: jacobi.cpp matrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
