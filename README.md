# Far Memory Latency Benchmark

This tool compares the performance of MemPageCache and MiraPageCache with simulated far memory latency. It helps demonstrate how Mira's cache eviction strategies perform under different far memory access latencies.

## Building

```bash
make clean
make latency
```

## Running

### Basic Run (3 iterations for each latency level)

```bash
./latency_benchmark
```

### Quick Run (2 iterations per latency level)

```bash
make latency_quick
```

### Thorough Run (10 iterations per latency level)

```bash
make latency_thorough
```

### Custom Run

```bash
./latency_benchmark [iterations] [page_size]

# Examples:
./latency_benchmark 5         # 5 iterations per latency level
./latency_benchmark 8 8192    # 8 iterations, 8192 bytes per page
```

## Latency Levels Tested

1. **Low Latency**: 100 microseconds - Simulates a high-performance CXL or RDMA connection
2. **Medium Latency**: 500 microseconds - Simulates typical networked memory or slower CXL
3. **High Latency**: 1 millisecond - Simulates far memory over a datacenter network

## Results

For each latency level, results are saved to separate files:
- `latency_benchmark_low.txt`
- `latency_benchmark_medium.txt`
- `latency_benchmark_high.txt`

Each file includes:
- Average, median, and 95th percentile times
- Performance difference percentages
- Hit/miss statistics for both cache implementations
- Eviction counts for MiraPageCache

## Interpreting Results

- **Negative difference**: MiraPageCache is faster
- **Positive difference**: MemPageCache is faster
- **Higher hit ratio**: More efficient cache utilization
- **Benefits of Mira should increase with latency**: As latency increases, the benefits of Mira's smart caching strategies should become more pronounced

## Tests Performed

1. **Insert Test**: Adding keys to the B+Tree
   - Tests how cache organization impacts write performance under latency

2. **Point Lookup Test**: Retrieving individual keys (with skewed access pattern)
   - 80% of lookups on 20% of data (hot data)
   - MiraPageCache optimizes for this pattern by prioritizing hot data

3. **Range Scan Test**: Retrieving ranges of consecutive keys
   - Tests sequential access optimization in high latency environments

4. **Mixed Workload Test**: Combination of lookups (70%), inserts (20%), and scans (10%)
   - Tests overall adaptability of each cache strategy under different latency conditions

## How This Relates to the Mira Paper

The Mira paper describes a far-memory system that uses program behavior to guide caching decisions. This benchmark demonstrates how such behavior-guided caching performs with different access latencies that might be encountered in real far-memory systems:

1. **Cache Section Separation**: Mira separates the cache into sections for different access patterns, which should particularly benefit higher latency scenarios

2. **Program-Guided Eviction**: As latency increases, making good eviction decisions becomes more critical, testing Mira's eviction hints

3. **Far Memory Simulation**: The paper describes performance with actual RDMA-based far memory - this benchmark uses simulated latency to approximate similar conditions

4. **Pattern-Specific Optimizations**: Tests whether optimizations for sequential, random, and stride patterns pay off more as latency increases