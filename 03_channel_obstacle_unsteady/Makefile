CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -pedantic

TARGET := channel_obstacle_unsteady
SRCS := main.cpp mesh.cpp boundary.cpp solver.cpp io.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) u.csv v.csv p.csv
	rm -rf output
