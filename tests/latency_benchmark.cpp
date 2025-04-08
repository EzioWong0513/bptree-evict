#include "../include/bptree/mira_page_cache_latency.h"
#include "../include/bptree/tree.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <random>
#include <chrono>
#include <thread>
#include <functional>
#include <numeric>
#include <algorithm>

// Structure to store statistics from multiple runs
struct BenchmarkStats {
    std::vector<double> mem_insert_times;
    std::vector<double> mira_insert_times;
    std::vector<double> mem_lookup_times;
    std::vector<double> mira_lookup_times;
    std::vector<double> mem_range_times;
    std::vector<double> mira_range_times;
    std::vector<double> mem_mixed_times;
    std::vector<double> mira_mixed_times;
    
    // Final statistics for caches
    uint64_t mem_hits = 0;
    uint64_t mem_misses = 0;
    uint64_t mira_hits = 0;
    uint64_t mira_misses = 0;
    uint64_t mira_evictions = 0;
};

// Function to measure execution time
template <typename Func>
double measure_time(Func func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    return duration.count();
}

// Calculate statistics from a vector of timing results
double calculate_average(const std::vector<double>& times) {
    if (times.empty()) return 0.0;
    return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
}

double calculate_median(std::vector<double> times) {
    if (times.empty()) return 0.0;
    std::sort(times.begin(), times.end());
    if (times.size() % 2 == 0) {
        return (times[times.size()/2 - 1] + times[times.size()/2]) / 2.0;
    } else {
        return times[times.size()/2];
    }
}

double calculate_percentile(std::vector<double> times, double percentile) {
    if (times.empty()) return 0.0;
    std::sort(times.begin(), times.end());
    int idx = static_cast<int>(times.size() * percentile / 100.0);
    return times[idx];
}

double calculate_stddev(const std::vector<double>& times) {
    if (times.size() <= 1) return 0.0;
    double avg = calculate_average(times);
    double sum_squared_diff = 0.0;
    for (double time : times) {
        double diff = time - avg;
        sum_squared_diff += diff * diff;
    }
    return std::sqrt(sum_squared_diff / (times.size() - 1));
}

// Run a single benchmark iteration with simulated latency
void run_benchmark_iteration(
    bptree::MemLatencyPageCache* mem_cache, 
    bptree::MiraLatencyPageCache* mira_cache,
    BenchmarkStats& stats,
    int num_inserts,
    int num_lookups,
    int num_range_scans,
    int range_size,
    int num_mixed_ops,
    bool verbose = false) 
{
    // Create B+Trees with the caches
    bptree::BTree<256, int, int> mem_tree(mem_cache);
    bptree::BTree<256, int, int> mira_tree(mira_cache);
    
    // Create random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Insert Test
    if (verbose) {
        std::cout << "Inserting " << num_inserts << " key-value pairs...\n";
    }
    
    // MemLatencyPageCache insert test
    if (verbose) std::cout << "  Testing MemLatencyPageCache inserts... ";
    double mem_insert_time = measure_time([&]() {
        for (int i = 0; i < num_inserts; i++) {
            mem_tree.insert(i, i * 100);
        }
    });
    if (verbose) std::cout << "done in " << mem_insert_time << " ms\n";
    stats.mem_insert_times.push_back(mem_insert_time);
    
    // MiraLatencyPageCache insert test
    if (verbose) std::cout << "  Testing MiraLatencyPageCache inserts... ";
    double mira_insert_time = measure_time([&]() {
        for (int i = 0; i < num_inserts; i++) {
            mira_tree.insert(i, i * 100);
            
            // Mark sequential pattern after a chunk of insertions - reduced frequency
            if (i > 0 && i % 1000 == 0) {
                mira_cache->assign_section(i, bptree::AccessPattern::SEQUENTIAL);
            }
        }
    });
    if (verbose) std::cout << "done in " << mira_insert_time << " ms\n";
    stats.mira_insert_times.push_back(mira_insert_time);
    
    // Reset cache stats after insertions
    mem_cache->reset_stats();
    mira_cache->reset_stats();
    
    // Point Lookup Test with skewed distribution
    if (verbose) {
        std::cout << "Performing " << num_lookups << " point lookups with skewed distribution...\n";
    }
    
    // Create skewed distribution (80% of lookups on 20% of data)
    const int HOT_DATA_SIZE = num_inserts / 5;  // 20% of data
    std::uniform_int_distribution<> hot_dist(0, HOT_DATA_SIZE - 1);
    std::uniform_int_distribution<> cold_dist(HOT_DATA_SIZE, num_inserts - 1);
    std::bernoulli_distribution hot_selector(0.8);  // 80% chance for hot data
    
    // For storing results
    std::vector<int> mem_values;
    std::vector<int> mira_values;
    
    // Track keys that will be looked up to mark them as important in Mira
    std::vector<int> keys_to_lookup;
    for (int i = 0; i < 100; i++) {
        int key = hot_selector(gen) ? hot_dist(gen) : cold_dist(gen);
        keys_to_lookup.push_back(key);
    }
    
    // Mark these keys as important in Mira (assign to random pattern section)
    for (int key : keys_to_lookup) {
        mira_cache->assign_section(key, bptree::AccessPattern::RANDOM);
    }
    
    // MemLatencyPageCache lookup test
    if (verbose) std::cout << "  Testing MemLatencyPageCache lookups... ";
    double mem_lookup_time = measure_time([&]() {
        for (int i = 0; i < num_lookups; i++) {
            int key = hot_selector(gen) ? hot_dist(gen) : cold_dist(gen);
            mem_values.clear();
            mem_tree.get_value(key, mem_values);
        }
    });
    if (verbose) std::cout << "done in " << mem_lookup_time << " ms\n";
    stats.mem_lookup_times.push_back(mem_lookup_time);
    
    // MiraLatencyPageCache lookup test
    if (verbose) std::cout << "  Testing MiraLatencyPageCache lookups... ";
    double mira_lookup_time = measure_time([&]() {
        for (int i = 0; i < num_lookups; i++) {
            int key = hot_selector(gen) ? hot_dist(gen) : cold_dist(gen);
            mira_values.clear();
            mira_tree.get_value(key, mira_values);
        }
    });
    if (verbose) std::cout << "done in " << mira_lookup_time << " ms\n";
    stats.mira_lookup_times.push_back(mira_lookup_time);
    
    // Reset cache stats after lookups
    mem_cache->reset_stats();
    mira_cache->reset_stats();
    
    // Range Scan Test
    if (verbose) {
        std::cout << "Performing " << num_range_scans << " range scans (each with " << range_size << " elements)...\n";
    }
    
    std::uniform_int_distribution<> range_start_dist(0, num_inserts - range_size - 1);
    
    // MemLatencyPageCache range scan test
    if (verbose) std::cout << "  Testing MemLatencyPageCache range scans... ";
    double mem_range_time = measure_time([&]() {
        for (int i = 0; i < num_range_scans; i++) {
            int start_key = range_start_dist(gen);
            int count = 0;
            
            for (auto it = mem_tree.begin(start_key); it != mem_tree.end() && count < range_size; ++it) {
                count++;
            }
        }
    });
    if (verbose) std::cout << "done in " << mem_range_time << " ms\n";
    stats.mem_range_times.push_back(mem_range_time);
    
    // MiraLatencyPageCache range scan test
    if (verbose) std::cout << "  Testing MiraLatencyPageCache range scans... ";
    
    // Mark sequence of keys as sequential for range scan (with reduced frequency)
    for (int i = 0; i < num_range_scans; i++) {
        int start_key = range_start_dist(gen);
        // Only mark every 10th key to reduce overhead
        for (int j = 0; j < range_size; j += 10) {
            mira_cache->assign_section(start_key + j, bptree::AccessPattern::SEQUENTIAL);
        }
    }
    
    double mira_range_time = measure_time([&]() {
        for (int i = 0; i < num_range_scans; i++) {
            int start_key = range_start_dist(gen);
            int count = 0;
            
            for (auto it = mira_tree.begin(start_key); it != mira_tree.end() && count < range_size; ++it) {
                count++;
            }
        }
    });
    if (verbose) std::cout << "done in " << mira_range_time << " ms\n";
    stats.mira_range_times.push_back(mira_range_time);
    
    // Reset cache stats after range scans
    mem_cache->reset_stats();
    mira_cache->reset_stats();
    
    // Mixed Workload Test
    if (verbose) {
        std::cout << "Performing " << num_mixed_ops << " mixed operations (70% lookups, 20% inserts, 10% scans)...\n";
    }
    
    std::uniform_int_distribution<> mixed_op_dist(1, 100);
    std::uniform_int_distribution<> insert_key_dist(num_inserts, num_inserts + num_mixed_ops);
    std::uniform_int_distribution<> scan_len_dist(10, 100);
    
    // MemLatencyPageCache mixed workload test
    if (verbose) std::cout << "  Testing MemLatencyPageCache mixed workload... ";
    double mem_mixed_time = measure_time([&]() {
        for (int i = 0; i < num_mixed_ops; i++) {
            int op = mixed_op_dist(gen);
            
            if (op <= 70) {
                // 70% chance: point lookup with skewed distribution
                int key = hot_selector(gen) ? hot_dist(gen) : cold_dist(gen);
                mem_values.clear();
                mem_tree.get_value(key, mem_values);
            } 
            else if (op <= 90) {
                // 20% chance: insert new key
                int key = insert_key_dist(gen);
                mem_tree.insert(key, key * 100);
            } 
            else {
                // 10% chance: small range scan
                int start_key = range_start_dist(gen);
                int scan_length = scan_len_dist(gen);
                int count = 0;
                
                for (auto it = mem_tree.begin(start_key); it != mem_tree.end() && count < scan_length; ++it) {
                    count++;
                }
            }
        }
    });
    if (verbose) std::cout << "done in " << mem_mixed_time << " ms\n";
    stats.mem_mixed_times.push_back(mem_mixed_time);
    
    // MiraLatencyPageCache mixed workload test
    if (verbose) std::cout << "  Testing MiraLatencyPageCache mixed workload... ";
    double mira_mixed_time = measure_time([&]() {
        for (int i = 0; i < num_mixed_ops; i++) {
            int op = mixed_op_dist(gen);
            
            if (op <= 70) {
                // 70% chance: point lookup with skewed distribution
                int key = hot_selector(gen) ? hot_dist(gen) : cold_dist(gen);
                mira_values.clear();
                mira_tree.get_value(key, mira_values);
                
                // Keep track of frequently accessed keys and prioritize them
                // Reduced frequency (every 50th operation instead of every 20th)
                if (i % 50 == 0) {
                    mira_cache->assign_section(key, bptree::AccessPattern::RANDOM);
                }
            } 
            else if (op <= 90) {
                // 20% chance: insert new key
                int key = insert_key_dist(gen);
                mira_tree.insert(key, key * 100);
                
                // Mark new key with appropriate pattern (reduced frequency)
                if (i % 10 == 0) {
                    if (i % 20 == 0) {
                        mira_cache->assign_section(key, bptree::AccessPattern::SEQUENTIAL);
                    } else {
                        mira_cache->assign_section(key, bptree::AccessPattern::STRIDE);
                    }
                }
            } 
            else {
                // 10% chance: small range scan
                int start_key = range_start_dist(gen);
                int scan_length = scan_len_dist(gen);
                int count = 0;
                
                // Mark range as sequential before scanning (reduced frequency)
                for (int j = 0; j < scan_length; j += 20) {
                    mira_cache->assign_section(start_key + j, bptree::AccessPattern::SEQUENTIAL);
                }
                
                for (auto it = mira_tree.begin(start_key); it != mira_tree.end() && count < scan_length; ++it) {
                    count++;
                }
            }
        }
    });
    if (verbose) std::cout << "done in " << mira_mixed_time << " ms\n";
    stats.mira_mixed_times.push_back(mira_mixed_time);
    
    // Capture final cache stats
    stats.mem_hits = mem_cache->get_hits();
    stats.mem_misses = mem_cache->get_misses();
    stats.mira_hits = mira_cache->get_hits();
    stats.mira_misses = mira_cache->get_misses();
    stats.mira_evictions = mira_cache->get_evictions();
}

// Run multiple benchmark iterations and generate a report for each latency level
void run_multi_benchmark(
    int num_iterations, 
    size_t page_size,
    bool warm_up = true) 
{
    std::cout << "====== Far Memory Latency Benchmark ======\n";
    std::cout << "Comparing MemLatencyPageCache vs MiraLatencyPageCache\n";
    std::cout << "Running " << num_iterations << " iterations per latency level\n";
    std::cout << "Page Size: " << page_size << " bytes\n";
    
    // Parameters for each test
    const int NUM_INSERTS = 100000;    
    const int NUM_LOOKUPS = 20000;
    const int NUM_RANGE_SCANS = 200;
    const int RANGE_SIZE = 500;
    const int NUM_MIXED_OPS = 20000;
    
    // Test with different latency levels
    std::vector<bptree::LatencyLevel> latency_levels = {
        bptree::LatencyLevel::LOW,      // 100 microseconds
        bptree::LatencyLevel::MEDIUM,   // 500 microseconds
        bptree::LatencyLevel::HIGH      // 1 millisecond
    };
    
    for (auto latency : latency_levels) {
        std::string latency_name;
        switch (latency) {
            case bptree::LatencyLevel::LOW: latency_name = "low"; break;
            case bptree::LatencyLevel::MEDIUM: latency_name = "medium"; break;
            case bptree::LatencyLevel::HIGH: latency_name = "high"; break;
            default: latency_name = "unknown";
        }
        
        std::cout << "\n\n===== Testing with " << latency_name << " latency =====\n";
        
        // Create output file for results
        std::string result_filename = "./latency_benchmark_" + latency_name + ".txt";
        std::ofstream result_file(result_filename);
        result_file << "Far Memory Latency Benchmark Results (" << latency_name << " latency)\n";
        result_file << "==============================================\n\n";
        result_file << "Configuration:\n";
        result_file << "  - Benchmark iterations: " << num_iterations << "\n";
        result_file << "  - Page size: " << page_size << " bytes\n";
        result_file << "  - Latency level: " << latency_name << "\n\n";
        
        // Storage for all benchmark results
        BenchmarkStats stats;
        
        // Run a warm-up iteration first if requested (results discarded)
        if (warm_up) {
            std::cout << "Running warm-up iteration... ";
            
            auto mem_cache = std::make_unique<bptree::MemLatencyPageCache>(page_size, latency);
            auto mira_cache = std::make_unique<bptree::MiraLatencyPageCache>(page_size, latency);
                
            BenchmarkStats warmup_stats;
            run_benchmark_iteration(
                mem_cache.get(), mira_cache.get(), warmup_stats,
                NUM_INSERTS, NUM_LOOKUPS, NUM_RANGE_SCANS, RANGE_SIZE, NUM_MIXED_OPS, false);
                
            std::cout << "done\n";
        }
        
        // Run the actual benchmark iterations
        for (int i = 0; i < num_iterations; i++) {
            std::cout << "Running iteration " << (i+1) << "/" << num_iterations << "... ";
            
            // Create fresh caches for each iteration
            auto mem_cache = std::make_unique<bptree::MemLatencyPageCache>(page_size, latency);
            auto mira_cache = std::make_unique<bptree::MiraLatencyPageCache>(page_size, latency);
                
            // Run the benchmark
            run_benchmark_iteration(
                mem_cache.get(), mira_cache.get(), stats,
                NUM_INSERTS, NUM_LOOKUPS, NUM_RANGE_SCANS, RANGE_SIZE, NUM_MIXED_OPS, false);
                
            std::cout << "done\n";
            
            // Sleep briefly between iterations to let system stabilize
            if (i < num_iterations - 1) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        }
        
        // Calculate statistics
        double mem_insert_avg = calculate_average(stats.mem_insert_times);
        double mira_insert_avg = calculate_average(stats.mira_insert_times);
        double mem_lookup_avg = calculate_average(stats.mem_lookup_times);
        double mira_lookup_avg = calculate_average(stats.mira_lookup_times);
        double mem_range_avg = calculate_average(stats.mem_range_times);
        double mira_range_avg = calculate_average(stats.mira_range_times);
        double mem_mixed_avg = calculate_average(stats.mem_mixed_times);
        double mira_mixed_avg = calculate_average(stats.mira_mixed_times);
        
        double mem_insert_median = calculate_median(stats.mem_insert_times);
        double mira_insert_median = calculate_median(stats.mira_insert_times);
        double mem_lookup_median = calculate_median(stats.mem_lookup_times);
        double mira_lookup_median = calculate_median(stats.mira_lookup_times);
        double mem_range_median = calculate_median(stats.mem_range_times);
        double mira_range_median = calculate_median(stats.mira_range_times);
        double mem_mixed_median = calculate_median(stats.mem_mixed_times);
        double mira_mixed_median = calculate_median(stats.mira_mixed_times);
        
        double mem_insert_p95 = calculate_percentile(stats.mem_insert_times, 95);
        double mira_insert_p95 = calculate_percentile(stats.mira_insert_times, 95);
        double mem_lookup_p95 = calculate_percentile(stats.mem_lookup_times, 95);
        double mira_lookup_p95 = calculate_percentile(stats.mira_lookup_times, 95);
        double mem_range_p95 = calculate_percentile(stats.mem_range_times, 95);
        double mira_range_p95 = calculate_percentile(stats.mira_range_times, 95);
        double mem_mixed_p95 = calculate_percentile(stats.mem_mixed_times, 95);
        double mira_mixed_p95 = calculate_percentile(stats.mira_mixed_times, 95);
        
        double mem_insert_stddev = calculate_stddev(stats.mem_insert_times);
        double mira_insert_stddev = calculate_stddev(stats.mira_insert_times);
        double mem_lookup_stddev = calculate_stddev(stats.mem_lookup_times);
        double mira_lookup_stddev = calculate_stddev(stats.mira_lookup_times);
        double mem_range_stddev = calculate_stddev(stats.mem_range_times);
        double mira_range_stddev = calculate_stddev(stats.mira_range_times);
        double mem_mixed_stddev = calculate_stddev(stats.mem_mixed_times);
        double mira_mixed_stddev = calculate_stddev(stats.mira_mixed_times);
        
        // Calculate percentage differences using median values (more robust than average)
        double insert_diff_pct = ((mira_insert_median - mem_insert_median) / mem_insert_median) * 100;
        double lookup_diff_pct = ((mira_lookup_median - mem_lookup_median) / mem_lookup_median) * 100;
        double range_diff_pct = ((mira_range_median - mem_range_median) / mem_range_median) * 100;
        double mixed_diff_pct = ((mira_mixed_median - mem_mixed_median) / mem_mixed_median) * 100;
        
        // Write results to file
        result_file << "\nPerformance Results (averaged over " << num_iterations << " iterations):\n";
        result_file << "--------------------------------------------------------\n";
        
        // Create a nice formatted comparison table
        result_file << std::fixed << std::setprecision(2);
        result_file << std::left << std::setw(20) << "Operation" 
                   << std::right << std::setw(12) << "MemCache" 
                   << std::setw(12) << "MiraCache" 
                   << std::setw(15) << "Diff (%)" 
                   << std::setw(11) << "MC StdDev" 
                   << std::setw(11) << "MC StdDev" << "\n";
        result_file << std::string(81, '-') << "\n";
        
        // Average results
        result_file << "AVERAGE TIMES (ms)\n";
        
        // Insert performance
        result_file << std::left << std::setw(20) << "Insert" 
                   << std::right << std::setw(12) << mem_insert_avg
                   << std::setw(12) << mira_insert_avg
                   << std::setw(15) << ((mira_insert_avg - mem_insert_avg) / mem_insert_avg) * 100
                   << std::setw(11) << mem_insert_stddev
                   << std::setw(11) << mira_insert_stddev << "\n";
        
        // Point lookup performance
        result_file << std::left << std::setw(20) << "Point Lookup" 
                   << std::right << std::setw(12) << mem_lookup_avg
                   << std::setw(12) << mira_lookup_avg
                   << std::setw(15) << ((mira_lookup_avg - mem_lookup_avg) / mem_lookup_avg) * 100
                   << std::setw(11) << mem_lookup_stddev
                   << std::setw(11) << mira_lookup_stddev << "\n";
        
        // Range scan performance
        result_file << std::left << std::setw(20) << "Range Scan" 
                   << std::right << std::setw(12) << mem_range_avg
                   << std::setw(12) << mira_range_avg
                   << std::setw(15) << ((mira_range_avg - mem_range_avg) / mem_range_avg) * 100
                   << std::setw(11) << mem_range_stddev
                   << std::setw(11) << mira_range_stddev << "\n";
        
        // Mixed workload performance
        result_file << std::left << std::setw(20) << "Mixed Workload" 
                   << std::right << std::setw(12) << mem_mixed_avg
                   << std::setw(12) << mira_mixed_avg
                   << std::setw(15) << ((mira_mixed_avg - mem_mixed_avg) / mem_mixed_avg) * 100
                   << std::setw(11) << mem_mixed_stddev
                   << std::setw(11) << mira_mixed_stddev << "\n\n";
        
        // Median results
        result_file << "MEDIAN TIMES (ms)\n";
        
        // Insert performance
        result_file << std::left << std::setw(20) << "Insert" 
                   << std::right << std::setw(12) << mem_insert_median
                   << std::setw(12) << mira_insert_median
                   << std::setw(15) << insert_diff_pct << "\n";
        
        // Point lookup performance
        result_file << std::left << std::setw(20) << "Point Lookup" 
                   << std::right << std::setw(12) << mem_lookup_median
                   << std::setw(12) << mira_lookup_median
                   << std::setw(15) << lookup_diff_pct << "\n";
        
        // Range scan performance
        result_file << std::left << std::setw(20) << "Range Scan" 
                   << std::right << std::setw(12) << mem_range_median
                   << std::setw(12) << mira_range_median
                   << std::setw(15) << range_diff_pct << "\n";
        
        // Mixed workload performance
        result_file << std::left << std::setw(20) << "Mixed Workload" 
                   << std::right << std::setw(12) << mem_mixed_median
                   << std::setw(12) << mira_mixed_median
                   << std::setw(15) << mixed_diff_pct << "\n\n";
        
        // Calculate cache hit ratios
        double mem_hit_ratio = (stats.mem_hits + stats.mem_misses > 0) ? 
            (stats.mem_hits * 100.0) / (stats.mem_hits + stats.mem_misses) : 0;
        double mira_hit_ratio = (stats.mira_hits + stats.mira_misses > 0) ?
            (stats.mira_hits * 100.0) / (stats.mira_hits + stats.mira_misses) : 0;
        
        // Cache statistics from the last iteration
        result_file << "Cache Statistics (from last iteration):\n";
        result_file << "-----------------------------------------\n";
        result_file << "MemLatencyPageCache:\n";
        result_file << "  - Cache hits: " << stats.mem_hits << "\n";
        result_file << "  - Cache misses: " << stats.mem_misses << "\n";
        result_file << "  - Hit ratio: " << mem_hit_ratio << "%\n\n";
        
        result_file << "MiraLatencyPageCache:\n";
        result_file << "  - Cache hits: " << stats.mira_hits << "\n";
        result_file << "  - Cache misses: " << stats.mira_misses << "\n";
        result_file << "  - Hit ratio: " << mira_hit_ratio << "%\n";
        result_file << "  - Evictions: " << stats.mira_evictions << "\n\n";
        
        // Print conclusion based on median values
        result_file << "Conclusion:\n";
        result_file << "-----------\n";
        if (mira_lookup_median < mem_lookup_median && mira_mixed_median < mem_mixed_median) {
            result_file << "With " << latency_name << " latency: MiraLatencyPageCache outperforms MemLatencyPageCache "
                       << "for both point lookups and mixed workloads, achieving " << std::abs(lookup_diff_pct) 
                       << "% faster lookups and " << std::abs(mixed_diff_pct) << "% faster mixed operations.\n";
        } else if (mira_lookup_median < mem_lookup_median) {
            result_file << "With " << latency_name << " latency: MiraLatencyPageCache outperforms MemLatencyPageCache "
                       << "for point lookups by " << std::abs(lookup_diff_pct) 
                       << "%, but is slower for other operations.\n";
        } else if (mem_lookup_median < mira_lookup_median) {
            result_file << "With " << latency_name << " latency: MemLatencyPageCache outperforms MiraLatencyPageCache "
                       << "for point lookups by " << std::abs(lookup_diff_pct) << "%.\n";
        }
        
        // If the hit ratio is significantly different, highlight that
        if (std::abs(mira_hit_ratio - mem_hit_ratio) > 5.0) {
            if (mira_hit_ratio > mem_hit_ratio) {
                result_file << "MiraLatencyPageCache achieves a " << (mira_hit_ratio - mem_hit_ratio) 
                           << "% higher hit ratio compared to MemLatencyPageCache.\n";
            } else {
                result_file << "MemLatencyPageCache achieves a " << (mem_hit_ratio - mira_hit_ratio) 
                           << "% higher hit ratio compared to MiraLatencyPageCache.\n";
            }
        }
        
        result_file.close();
        
        // Print summary to console
        std::cout << "\n========== " << latency_name << " Latency Summary ==========\n";
        std::cout << std::fixed << std::setprecision(2);
        
        std::cout << "MemLatencyPageCache vs MiraLatencyPageCache performance comparison (median values):\n";
        std::cout << "  - Insert: " << mem_insert_median << " ms vs " << mira_insert_median 
                 << " ms (" << (insert_diff_pct > 0 ? "+" : "") << insert_diff_pct << "%)\n";
        
        std::cout << "  - Point Lookup: " << mem_lookup_median << " ms vs " << mira_lookup_median 
                 << " ms (" << (lookup_diff_pct > 0 ? "+" : "") << lookup_diff_pct << "%)\n";
        
        std::cout << "  - Range Scan: " << mem_range_median << " ms vs " << mira_range_median 
                 << " ms (" << (range_diff_pct > 0 ? "+" : "") << range_diff_pct << "%)\n";
        
        std::cout << "  - Mixed Workload: " << mem_mixed_median << " ms vs " << mira_mixed_median 
                 << " ms (" << (mixed_diff_pct > 0 ? "+" : "") << mixed_diff_pct << "%)\n";
        
        // Report cache hit ratios
        std::cout << "\nMemLatencyPageCache hit ratio: " << mem_hit_ratio << "%\n";
        std::cout << "MiraLatencyPageCache hit ratio: " << mira_hit_ratio << "%\n";
        
        std::cout << "\nDetailed results written to: " << result_filename << "\n";
        std::cout << "=====================================\n";
    }
}

int main(int argc, char* argv[]) {
    try {
        // Default values
        int num_iterations = 3;
        size_t page_size = 4096;
        bool warm_up = true;
        
        // Parse command line arguments if provided
        if (argc > 1) {
            num_iterations = std::stoi(argv[1]);
        }
        
        if (argc > 2) {
            page_size = std::stoul(argv[2]);
        }
        
        // Run benchmark with specified iterations
        run_multi_benchmark(num_iterations, page_size, warm_up);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}