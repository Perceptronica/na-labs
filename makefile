CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -O2
TARGET = calc
SRCS = main.cpp src/linalg.cpp src/matrix.cpp src/wrappers.cpp
OBJS = $(SRCS:.cpp=.o)
HEADERS = lib/linalg.h lib/matrix.h lib/wrappers.h
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
