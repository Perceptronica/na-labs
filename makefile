CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2
TARGET = matrix_program
SRCS = main.cpp matrix.cpp
OBJS = $(SRCS:.cpp=.o)
HEADERS = matrix.h
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
