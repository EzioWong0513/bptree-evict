CXX = g++

CXXFLAGS = -std=c++17 -I../include -I/usr/include

LIBS = -lpthread -lboost_thread -lboost_system

TARGET = main
MIRA_TARGET = mira_test
BENCHMARK_TARGET = mira_compare_benchmark
LATENCY_TARGET = latency_benchmark
LIMITED_LATENCY_TARGET = limited_latency_benchmark

SRCS = tests/main.cpp src/heap_page_cache.cpp src/heap_file.cpp
MIRA_SRCS = tests/mira_test.cpp src/heap_page_cache.cpp src/heap_file.cpp
BENCHMARK_SRCS = tests/mira_compare_benchmark.cpp src/heap_page_cache.cpp src/heap_file.cpp
LATENCY_SRCS = tests/latency_benchmark.cpp src/heap_page_cache.cpp src/heap_file.cpp
LIMITED_LATENCY_SRCS = tests/limited_latency_benchmark.cpp src/heap_page_cache.cpp src/heap_file.cpp

OBJS = $(SRCS:.cpp=.o)
MIRA_OBJS = $(MIRA_SRCS:.cpp=.o)
BENCHMARK_OBJS = $(BENCHMARK_SRCS:.cpp=.o)
LATENCY_OBJS = $(LATENCY_SRCS:.cpp=.o)
LIMITED_LATENCY_OBJS = $(LIMITED_LATENCY_SRCS:.cpp=.o)

all: $(TARGET) $(MIRA_TARGET) $(BENCHMARK_TARGET) $(LATENCY_TARGET) $(LIMITED_LATENCY_TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS) $(LIBS)

$(MIRA_TARGET): $(MIRA_OBJS)
	$(CXX) $(CXXFLAGS) -o $(MIRA_TARGET) $(MIRA_SRCS) $(LIBS)

$(BENCHMARK_TARGET): $(BENCHMARK_OBJS)
	$(CXX) $(CXXFLAGS) -o $(BENCHMARK_TARGET) $(BENCHMARK_SRCS) $(LIBS)

$(LATENCY_TARGET): $(LATENCY_OBJS)
	$(CXX) $(CXXFLAGS) -o $(LATENCY_TARGET) $(LATENCY_SRCS) $(LIBS)

$(LIMITED_LATENCY_TARGET): $(LIMITED_LATENCY_OBJS)
	$(CXX) $(CXXFLAGS) -o $(LIMITED_LATENCY_TARGET) $(LIMITED_LATENCY_SRCS) $(LIBS)

tests/%.o: tests/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

tests: all
	./$(TARGET)

mira_tests: $(MIRA_TARGET)
	./$(MIRA_TARGET)

benchmark: $(BENCHMARK_TARGET)
	./$(BENCHMARK_TARGET)

latency: $(LATENCY_TARGET)
	./$(LATENCY_TARGET)

limited_latency: $(LIMITED_LATENCY_TARGET)
	./$(LIMITED_LATENCY_TARGET)

# Run latency tests with specific parameters
latency_quick: $(LATENCY_TARGET)
	./$(LATENCY_TARGET) 2

latency_thorough: $(LATENCY_TARGET)
	./$(LATENCY_TARGET) 10

# Run limited latency tests with specific parameters
limited_latency_quick: $(LIMITED_LATENCY_TARGET)
	./$(LIMITED_LATENCY_TARGET) 2

limited_latency_thorough: $(LIMITED_LATENCY_TARGET)
	./$(LIMITED_LATENCY_TARGET) 10

clean:
	rm -f $(TARGET) $(MIRA_TARGET) $(BENCHMARK_TARGET) $(LATENCY_TARGET) $(LIMITED_LATENCY_TARGET) $(OBJS) $(MIRA_OBJS) $(BENCHMARK_OBJS) $(LATENCY_OBJS) $(LIMITED_LATENCY_OBJS)
	rm -rf tmp/*
	rm -f profile.pdf profile.svg profile.prof
	rm -f *.txt

.PHONY: all clean tests mira_tests benchmark latency limited_latency latency_quick latency_thorough limited_latency_quick limited_latency_thorough