CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = channel_flow
SRCS   = main.cpp solver.cpp boundary.cpp io.cpp
OBJS   = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET) result.csv result.vtk

.PHONY: clean
