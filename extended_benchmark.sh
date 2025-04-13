#!/bin/bash
# Script to build and run the extended benchmark

# Set the number of iterations 
if [ -z "$1" ]; then
    ITERATIONS=50
else
    ITERATIONS=$1
fi

echo "===== Building Extended Benchmark ====="
# Create build directory
mkdir -p build
cd build

# Compile the enhanced benchmark
g++ -o extended_benchmark ../enhanced_benchmark.cpp ../src/heap_page_cache.cpp ../src/heap_file.cpp -I../include -std=c++17 -lpthread -lboost_thread -lboost_system

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "===== Running Extended Benchmark with $ITERATIONS iterations ====="
# Run the benchmark
./extended_benchmark $ITERATIONS

if [ $? -ne 0 ]; then
    echo "Benchmark execution failed!"
    exit 1
fi

cd ..

# Copy CSV files to the main directory
cp build/*.csv ..

echo "===== Generating Plots ====="
# Check if Python is available
if command -v python3 &>/dev/null; then
    # Check if required Python packages are installed
    python3 -c "import pandas, matplotlib, seaborn" 2>/dev/null
    if [ $? -eq 0 ]; then
        # Run the plotting script
        python3 plot_benchmark_results.py
    else
        echo "Required Python packages not found. Please install them with:"
        echo "pip install pandas matplotlib seaborn"
        echo "Then run: python3 plot_benchmark_results.py"
    fi
else
    echo "Python 3 not found. Please install Python and run:"
    echo "python3 plot_benchmark_results.py"
fi

echo "===== Benchmark Complete ====="
echo "Results can be found in the benchmark_results directory"
echo "Plots can be found in the benchmark_plots directory"