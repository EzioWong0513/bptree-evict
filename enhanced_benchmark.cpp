#include "../include/bptree/mira_page_cache_latency.h"
#include "../include/bptree/tree.h"
#include "../include/bptree/limited_mem_latency_page_cache.h"
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
#include <map>
#include <cmath>
#include <filesystem>

// Structure to store detailed statistics from multiple runs
struct DetailedBenchmarkStats {
    // Raw timing data for each operation and each iteration
    std::vector<double> mem_insert_times;
    std::vector<double> mira_insert_times;
    std::vector<double> mem_lookup_times;
    std::vector<double> mira_lookup_times;
    std::vector<double> mem_range_times;
    std::vector<double> mira_range_times;
    std::vector<double> mem_mixed_times;
    std::vector<double> mira_mixed_times;
    
    // Cache performance metrics for each iteration
    std::vector<uint64_t> mem_hits_per_iter;
    std::vector<uint64_t> mem_misses_per_iter;
    std::vector<uint64_t> mem_evictions_per_iter;
    std::vector<uint64_t> mira_hits_per_iter;
    std::vector<uint64_t> mira_misses_per_iter;
    std::vector<uint64_t> mira_evictions_per_iter;
    
    // Final statistics for caches (last iteration)
    uint64_t mem_hits = 0;
    uint64_t mem_misses = 0;
    uint64_t mem_evictions = 0;
    uint64_t mira_hits = 0;
    uint64_t mira_misses = 0;
    uint64_t mira_evictions = 0;
};

// Map to store stats for each latency level
std::map<bptree::LatencyLevel, DetailedBenchmarkStats> all_stats;

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
    if (idx >= times.size()) idx = times.size() - 1;
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
    bptree::LimitedMemLatencyPageCache* mem_cache,
    bptree::MiraLatencyPageCache* mira_cache,
    DetailedBenchmarkStats& stats,
    int num_inserts,
    int num_lookups,
    int num_range_scans,
    int range_size,
    int num_mixed_ops,
    int iteration,
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
    
    // Reset cache stats before test
    mem_cache->reset_stats();
    mira_cache->reset_stats();
    
    // MemLatencyPageCache insert test
    if (verbose) std::cout << "  Testing MemLatencyPageCache inserts... ";
    double mem_insert_time = measure_time([&]() {
        for (int i = 0; i < num_inserts; i++) {
            mem_tree.insert(i, i * 100);
        }
    });
    if (verbose) std::cout << "done in " << mem_insert_time << " ms\n";
    stats.mem_insert_times.push_back(mem_insert_time);
    
    // Save cache stats after mem cache test
    stats.mem_hits_per_iter.push_back(mem_cache->get_hits());
    stats.mem_misses_per_iter.push_back(mem_cache->get_misses());
    stats.mem_evictions_per_iter.push_back(mem_cache->get_evictions());
    
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
    
    // Save cache stats after mira cache test
    stats.mira_hits_per_iter.push_back(mira_cache->get_hits());
    stats.mira_misses_per_iter.push_back(mira_cache->get_misses());
    stats.mira_evictions_per_iter.push_back(mira_cache->get_evictions());
    
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
    
    // Update cache stats after mem lookup test
    stats.mem_hits_per_iter.push_back(mem_cache->get_hits());
    stats.mem_misses_per_iter.push_back(mem_cache->get_misses());
    stats.mem_evictions_per_iter.push_back(mem_cache->get_evictions());
    
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
    
    // Update cache stats after mira lookup test
    stats.mira_hits_per_iter.push_back(mira_cache->get_hits());
    stats.mira_misses_per_iter.push_back(mira_cache->get_misses());
    stats.mira_evictions_per_iter.push_back(mira_cache->get_evictions());
    
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
    
    // Update cache stats after mem range scan test
    stats.mem_hits_per_iter.push_back(mem_cache->get_hits());
    stats.mem_misses_per_iter.push_back(mem_cache->get_misses());
    stats.mem_evictions_per_iter.push_back(mem_cache->get_evictions());
    
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
    
    // Update cache stats after mira range scan test
    stats.mira_hits_per_iter.push_back(mira_cache->get_hits());
    stats.mira_misses_per_iter.push_back(mira_cache->get_misses());
    stats.mira_evictions_per_iter.push_back(mira_cache->get_evictions());
    
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
    
    // Update cache stats after mem mixed workload test
    stats.mem_hits_per_iter.push_back(mem_cache->get_hits());
    stats.mem_misses_per_iter.push_back(mem_cache->get_misses());
    stats.mem_evictions_per_iter.push_back(mem_cache->get_evictions());
    
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
    
    // Update cache stats after mira mixed workload test
    stats.mira_hits_per_iter.push_back(mira_cache->get_hits());
    stats.mira_misses_per_iter.push_back(mira_cache->get_misses());
    stats.mira_evictions_per_iter.push_back(mira_cache->get_evictions());
    
    // Capture final cache stats (only on the last iteration)
    stats.mem_hits = mem_cache->get_hits();
    stats.mem_misses = mem_cache->get_misses();
    stats.mem_evictions = mem_cache->get_evictions();
    stats.mira_hits = mira_cache->get_hits();
    stats.mira_misses = mira_cache->get_misses();
    stats.mira_evictions = mira_cache->get_evictions();
    
    // Output progress for long-running benchmarks
    if (iteration % 5 == 0) {
        std::cout << "  Completed iteration " << iteration << std::endl;
    }
}

// Generate CSV files for plotting
void generate_csv_files(const std::map<bptree::LatencyLevel, DetailedBenchmarkStats>& all_stats) {
    // Map latency level to name
    std::map<bptree::LatencyLevel, std::string> latency_names = {
        {bptree::LatencyLevel::LOW, "low"},
        {bptree::LatencyLevel::MEDIUM, "medium"},
        {bptree::LatencyLevel::HIGH, "high"}
    };
    
    // Generate CSV for lookup times
    std::ofstream lookup_csv("lookup_times.csv");
    lookup_csv << "Iteration,LatencyLevel,MemCache,MiraCache,Difference(%)\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        for (size_t i = 0; i < stats.mem_lookup_times.size(); ++i) {
            double mem_time = stats.mem_lookup_times[i];
            double mira_time = stats.mira_lookup_times[i];
            double diff_pct = ((mira_time - mem_time) / mem_time) * 100;
            
            lookup_csv << (i+1) << "," << latency_name << "," 
                      << mem_time << "," << mira_time << "," 
                      << diff_pct << "\n";
        }
    }
    lookup_csv.close();
    
    // Generate CSV for range scan times
    std::ofstream range_csv("range_scan_times.csv");
    range_csv << "Iteration,LatencyLevel,MemCache,MiraCache,Difference(%)\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        for (size_t i = 0; i < stats.mem_range_times.size(); ++i) {
            double mem_time = stats.mem_range_times[i];
            double mira_time = stats.mira_range_times[i];
            double diff_pct = ((mira_time - mem_time) / mem_time) * 100;
            
            range_csv << (i+1) << "," << latency_name << "," 
                     << mem_time << "," << mira_time << "," 
                     << diff_pct << "\n";
        }
    }
    range_csv.close();
    
    // Generate CSV for mixed workload times
    std::ofstream mixed_csv("mixed_workload_times.csv");
    mixed_csv << "Iteration,LatencyLevel,MemCache,MiraCache,Difference(%)\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        for (size_t i = 0; i < stats.mem_mixed_times.size(); ++i) {
            double mem_time = stats.mem_mixed_times[i];
            double mira_time = stats.mira_mixed_times[i];
            double diff_pct = ((mira_time - mem_time) / mem_time) * 100;
            
            mixed_csv << (i+1) << "," << latency_name << "," 
                     << mem_time << "," << mira_time << "," 
                     << diff_pct << "\n";
        }
    }
    mixed_csv.close();
    
    // Generate CSV for hit ratios
    std::ofstream hit_ratio_csv("hit_ratios.csv");
    hit_ratio_csv << "Operation,LatencyLevel,MemCacheHitRatio,MiraCacheHitRatio\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        
        // Insert operation hit ratios
        double mem_insert_hits = stats.mem_hits_per_iter[0];
        double mem_insert_total = mem_insert_hits + stats.mem_misses_per_iter[0];
        double mira_insert_hits = stats.mira_hits_per_iter[0];
        double mira_insert_total = mira_insert_hits + stats.mira_misses_per_iter[0];
        
        double mem_insert_ratio = mem_insert_total > 0 ? (mem_insert_hits / mem_insert_total) * 100 : 0;
        double mira_insert_ratio = mira_insert_total > 0 ? (mira_insert_hits / mira_insert_total) * 100 : 0;
        
        hit_ratio_csv << "Insert," << latency_name << "," 
                     << mem_insert_ratio << "," << mira_insert_ratio << "\n";
        
        // Lookup operation hit ratios
        double mem_lookup_hits = stats.mem_hits_per_iter[1];
        double mem_lookup_total = mem_lookup_hits + stats.mem_misses_per_iter[1];
        double mira_lookup_hits = stats.mira_hits_per_iter[1];
        double mira_lookup_total = mira_lookup_hits + stats.mira_misses_per_iter[1];
        
        double mem_lookup_ratio = mem_lookup_total > 0 ? (mem_lookup_hits / mem_lookup_total) * 100 : 0;
        double mira_lookup_ratio = mira_lookup_total > 0 ? (mira_lookup_hits / mira_lookup_total) * 100 : 0;
        
        hit_ratio_csv << "Lookup," << latency_name << "," 
                     << mem_lookup_ratio << "," << mira_lookup_ratio << "\n";
        
        // Range scan operation hit ratios
        double mem_range_hits = stats.mem_hits_per_iter[2];
        double mem_range_total = mem_range_hits + stats.mem_misses_per_iter[2];
        double mira_range_hits = stats.mira_hits_per_iter[2];
        double mira_range_total = mira_range_hits + stats.mira_misses_per_iter[2];
        
        double mem_range_ratio = mem_range_total > 0 ? (mem_range_hits / mem_range_total) * 100 : 0;
        double mira_range_ratio = mira_range_total > 0 ? (mira_range_hits / mira_range_total) * 100 : 0;
        
        hit_ratio_csv << "RangeScan," << latency_name << "," 
                     << mem_range_ratio << "," << mira_range_ratio << "\n";
        
        // Mixed workload operation hit ratios
        double mem_mixed_hits = stats.mem_hits_per_iter[3];
        double mem_mixed_total = mem_mixed_hits + stats.mem_misses_per_iter[3];
        double mira_mixed_hits = stats.mira_hits_per_iter[3];
        double mira_mixed_total = mira_mixed_hits + stats.mira_misses_per_iter[3];
        
        double mem_mixed_ratio = mem_mixed_total > 0 ? (mem_mixed_hits / mem_mixed_total) * 100 : 0;
        double mira_mixed_ratio = mira_mixed_total > 0 ? (mira_mixed_hits / mira_mixed_total) * 100 : 0;
        
        hit_ratio_csv << "Mixed," << latency_name << "," 
                     << mem_mixed_ratio << "," << mira_mixed_ratio << "\n";
    }
    hit_ratio_csv.close();
    
    // Generate CSV for summary statistics
    std::ofstream summary_csv("summary_stats.csv");
    summary_csv << "Operation,LatencyLevel,MemCacheAvg,MiraCacheAvg,MemCacheMedian,MiraCacheMedian,DifferenceAvg(%),DifferenceMedian(%)\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        
        // Insert operation summary
        double mem_insert_avg = calculate_average(stats.mem_insert_times);
        double mira_insert_avg = calculate_average(stats.mira_insert_times);
        double mem_insert_median = calculate_median(stats.mem_insert_times);
        double mira_insert_median = calculate_median(stats.mira_insert_times);
        double insert_avg_diff = ((mira_insert_avg - mem_insert_avg) / mem_insert_avg) * 100;
        double insert_median_diff = ((mira_insert_median - mem_insert_median) / mem_insert_median) * 100;
        
        summary_csv << "Insert," << latency_name << "," 
                   << mem_insert_avg << "," << mira_insert_avg << "," 
                   << mem_insert_median << "," << mira_insert_median << "," 
                   << insert_avg_diff << "," << insert_median_diff << "\n";
        
        // Lookup operation summary
        double mem_lookup_avg = calculate_average(stats.mem_lookup_times);
        double mira_lookup_avg = calculate_average(stats.mira_lookup_times);
        double mem_lookup_median = calculate_median(stats.mem_lookup_times);
        double mira_lookup_median = calculate_median(stats.mira_lookup_times);
        double lookup_avg_diff = ((mira_lookup_avg - mem_lookup_avg) / mem_lookup_avg) * 100;
        double lookup_median_diff = ((mira_lookup_median - mem_lookup_median) / mem_lookup_median) * 100;
        
        summary_csv << "Lookup," << latency_name << "," 
                   << mem_lookup_avg << "," << mira_lookup_avg << "," 
                   << mem_lookup_median << "," << mira_lookup_median << "," 
                   << lookup_avg_diff << "," << lookup_median_diff << "\n";
        
        // Range scan operation summary
        double mem_range_avg = calculate_average(stats.mem_range_times);
        double mira_range_avg = calculate_average(stats.mira_range_times);
        double mem_range_median = calculate_median(stats.mem_range_times);
        double mira_range_median = calculate_median(stats.mira_range_times);
        double range_avg_diff = ((mira_range_avg - mem_range_avg) / mem_range_avg) * 100;
        double range_median_diff = ((mira_range_median - mem_range_median) / mem_range_median) * 100;
        
        summary_csv << "RangeScan," << latency_name << "," 
                   << mem_range_avg << "," << mira_range_avg << "," 
                   << mem_range_median << "," << mira_range_median << "," 
                   << range_avg_diff << "," << range_median_diff << "\n";
        
        // Mixed workload operation summary
        double mem_mixed_avg = calculate_average(stats.mem_mixed_times);
        double mira_mixed_avg = calculate_average(stats.mira_mixed_times);
        double mem_mixed_median = calculate_median(stats.mem_mixed_times);
        double mira_mixed_median = calculate_median(stats.mira_mixed_times);
        double mixed_avg_diff = ((mira_mixed_avg - mem_mixed_avg) / mem_mixed_avg) * 100;
        double mixed_median_diff = ((mira_mixed_median - mem_mixed_median) / mem_mixed_median) * 100;
        
        summary_csv << "Mixed," << latency_name << "," 
                   << mem_mixed_avg << "," << mira_mixed_avg << "," 
                   << mem_mixed_median << "," << mira_mixed_median << "," 
                   << mixed_avg_diff << "," << mixed_median_diff << "\n";
    }
    summary_csv.close();
    
    // Generate CSV for eviction statistics
    std::ofstream eviction_csv("eviction_stats.csv");
    eviction_csv << "LatencyLevel,MemCacheEvictions,MiraCacheEvictions,EvictionRatio\n";
    
    for (const auto& [latency, stats] : all_stats) {
        std::string latency_name = latency_names[latency];
        double eviction_ratio = stats.mem_evictions > 0 ? 
            static_cast<double>(stats.mira_evictions) / stats.mem_evictions : 0;
        
        eviction_csv << latency_name << "," 
                    << stats.mem_evictions << "," 
                    << stats.mira_evictions << "," 
                    << eviction_ratio << "\n";
    }
    eviction_csv.close();
    
    std::cout << "CSV files generated for plotting:\n";
    std::cout << "  - lookup_times.csv\n";
    std::cout << "  - range_scan_times.csv\n";
    std::cout << "  - mixed_workload_times.csv\n";
    std::cout << "  - hit_ratios.csv\n";
    std::cout << "  - summary_stats.csv\n";
    std::cout << "  - eviction_stats.csv\n";
}

// Run multiple benchmark iterations and generate a report for each latency level
void run_multi_benchmark(
    int num_iterations, 
    size_t page_size,
    bool generate_plots = true,
    bool warm_up = true) 
{
    std::cout << "====== Far Memory Latency Benchmark (Extended Version) ======\n";
    std::cout << "Comparing LimitedMemLatencyPageCache vs MiraLatencyPageCache\n";
    std::cout << "Running " << num_iterations << " iterations per latency level\n";
    std::cout << "Page Size: " << page_size << " bytes\n";
    
    // Create results directory if it doesn't exist
    std::filesystem::create_directory("benchmark_results");
    
    // Parameters for each test - adjusted for better granularity
    const int NUM_INSERTS = 1000000;   
    const int NUM_LOOKUPS = 500000;   
    const int NUM_RANGE_SCANS = 200;
    const int RANGE_SIZE = 5000;
    const int NUM_MIXED_OPS = 20000;
    
    const size_t CACHE_SIZE_BYTES = 64 * 1024;  
    std::cout << "Using restricted cache size: " << (CACHE_SIZE_BYTES / 1024) << " KB\n";
    
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
        std::string result_filename = "benchmark_results/latency_benchmark_" + latency_name + "_extended.txt";
        std::ofstream result_file(result_filename);
        result_file << "Far Memory Latency Benchmark Extended Results (" << latency_name << " latency)\n";
        result_file << "==========================================================================\n\n";
        result_file << "Configuration:\n";
        result_file << "  - Benchmark iterations: " << num_iterations << "\n";
        result_file << "  - Page size: " << page_size << " bytes\n";
        result_file << "  - Cache size: " << (CACHE_SIZE_BYTES / 1024) << " KB\n";
        result_file << "  - Latency level: " << latency_name << "\n\n";
        
        // Storage for all benchmark results
        DetailedBenchmarkStats stats;
        
        // Run a warm-up iteration first if requested (results discarded)
        if (warm_up) {
            std::cout << "Running warm-up iteration... ";
            
            // Create caches with LIMITED CACHE SIZE for both
            auto mem_cache = std::make_unique<bptree::LimitedMemLatencyPageCache>(
                page_size, latency, CACHE_SIZE_BYTES);
            auto mira_cache = std::make_unique<bptree::MiraLatencyPageCache>(
                page_size, latency, CACHE_SIZE_BYTES);
                
            DetailedBenchmarkStats warmup_stats;
            run_benchmark_iteration(
                mem_cache.get(), mira_cache.get(), warmup_stats,
                NUM_INSERTS, NUM_LOOKUPS, NUM_RANGE_SCANS, RANGE_SIZE, NUM_MIXED_OPS, 0, false);
                
            std::cout << "done\n";
        }
        
        // Run the actual benchmark iterations
        std::cout << "Running " << num_iterations << " benchmark iterations with " << latency_name << " latency...\n";
        
        for (int i = 0; i < num_iterations; i++) {
            // Create fresh caches for each iteration WITH LIMITED CACHE SIZE for both
            auto mem_cache = std::make_unique<bptree::LimitedMemLatencyPageCache>(
                page_size, latency, CACHE_SIZE_BYTES);
            auto mira_cache = std::make_unique<bptree::MiraLatencyPageCache>(
                page_size, latency, CACHE_SIZE_BYTES);
                
            // Run the benchmark
            run_benchmark_iteration(
                mem_cache.get(), mira_cache.get(), stats,
                NUM_INSERTS, NUM_LOOKUPS, NUM_RANGE_SCANS, RANGE_SIZE, NUM_MIXED_OPS, i+1, false);
                
            // Sleep briefly between iterations to let system stabilize
            if (i < num_iterations - 1) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        
        std::cout << "Completed " << num_iterations << " iterations.\n";
        
        // Store stats in global map for later CSV generation
        all_stats[latency] = stats;
        
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
        
        double mem_lookup_p95 = calculate_percentile(stats.mem_lookup_times, 95);
        double mira_lookup_p95 = calculate_percentile(stats.mira_lookup_times, 95);
        double mem_range_p95 = calculate_percentile(stats.mem_range_times, 95);
        double mira_range_p95 = calculate_percentile(stats.mira_range_times, 95);
        
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
        
        // Calculate confidence intervals (95%)
        double lookup_ci_width_mem = 1.96 * mem_lookup_stddev / std::sqrt(stats.mem_lookup_times.size());
        double lookup_ci_width_mira = 1.96 * mira_lookup_stddev / std::sqrt(stats.mira_lookup_times.size());
        
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
        
        // 95th percentile results (for key operations)
        result_file << "95th PERCENTILE TIMES (ms)\n";
        result_file << std::left << std::setw(20) << "Point Lookup" 
                   << std::right << std::setw(12) << mem_lookup_p95
                   << std::setw(12) << mira_lookup_p95
                   << std::setw(15) << ((mira_lookup_p95 - mem_lookup_p95) / mem_lookup_p95) * 100 << "\n";
        
        result_file << std::left << std::setw(20) << "Range Scan" 
                   << std::right << std::setw(12) << mem_range_p95
                   << std::setw(12) << mira_range_p95
                   << std::setw(15) << ((mira_range_p95 - mem_range_p95) / mem_range_p95) * 100 << "\n\n";
        
        // Confidence interval for lookup performance (most critical operation)
        result_file << "CONFIDENCE INTERVALS (95%)\n";
        result_file << std::left << std::setw(20) << "Point Lookup" 
                   << std::right << std::setw(12) << mem_lookup_avg << " ± " << lookup_ci_width_mem
                   << std::setw(12) << mira_lookup_avg << " ± " << lookup_ci_width_mira << "\n\n";
        
        // Calculate cache hit ratios
        double mem_hit_ratio = (stats.mem_hits + stats.mem_misses > 0) ? 
            (stats.mem_hits * 100.0) / (stats.mem_hits + stats.mem_misses) : 0;
        double mira_hit_ratio = (stats.mira_hits + stats.mira_misses > 0) ?
            (stats.mira_hits * 100.0) / (stats.mira_hits + stats.mira_misses) : 0;
        
        // Cache statistics from the last iteration
        result_file << "Cache Statistics (from last iteration):\n";
        result_file << "-----------------------------------------\n";
        result_file << "LimitedMemLatencyPageCache:\n";
        result_file << "  - Cache hits: " << stats.mem_hits << "\n";
        result_file << "  - Cache misses: " << stats.mem_misses << "\n";
        result_file << "  - Hit ratio: " << mem_hit_ratio << "%\n";
        result_file << "  - Evictions: " << stats.mem_evictions << "\n\n";
        
        result_file << "MiraLatencyPageCache:\n";
        result_file << "  - Cache hits: " << stats.mira_hits << "\n";
        result_file << "  - Cache misses: " << stats.mira_misses << "\n";
        result_file << "  - Hit ratio: " << mira_hit_ratio << "%\n";
        result_file << "  - Evictions: " << stats.mira_evictions << "\n\n";
        
        // Statistical significance assessment
        bool is_significant = std::abs(lookup_diff_pct) > (2 * lookup_ci_width_mem / mem_lookup_avg * 100);
        
        // Print conclusion based on median values and statistical significance
        result_file << "Conclusion:\n";
        result_file << "-----------\n";
        
        if (is_significant) {
            if (mira_lookup_median < mem_lookup_median) {
                result_file << "With " << latency_name << " latency: MiraLatencyPageCache significantly outperforms "
                           << "LimitedMemLatencyPageCache for point lookups by " << std::abs(lookup_diff_pct) 
                           << "% (statistically significant).\n";
            } else if (mem_lookup_median < mira_lookup_median) {
                result_file << "With " << latency_name << " latency: LimitedMemLatencyPageCache significantly outperforms "
                           << "MiraLatencyPageCache for point lookups by " << std::abs(lookup_diff_pct) 
                           << "% (statistically significant).\n";
            }
        } else {
            if (mira_lookup_median < mem_lookup_median) {
                result_file << "With " << latency_name << " latency: MiraLatencyPageCache appears to outperform "
                           << "LimitedMemLatencyPageCache for point lookups by " << std::abs(lookup_diff_pct) 
                           << "%, but the difference is not statistically significant.\n";
            } else if (mem_lookup_median < mira_lookup_median) {
                result_file << "With " << latency_name << " latency: LimitedMemLatencyPageCache appears to outperform "
                           << "MiraLatencyPageCache for point lookups by " << std::abs(lookup_diff_pct) 
                           << "%, but the difference is not statistically significant.\n";
            } else {
                result_file << "With " << latency_name << " latency: There is no significant performance difference "
                           << "between LimitedMemLatencyPageCache and MiraLatencyPageCache for point lookups.\n";
            }
        }
        
        // Compare cache hit ratios
        if (std::abs(mira_hit_ratio - mem_hit_ratio) > 1.0) {
            if (mira_hit_ratio > mem_hit_ratio) {
                result_file << "MiraLatencyPageCache achieves a " << (mira_hit_ratio - mem_hit_ratio) 
                           << "% higher hit ratio compared to LimitedMemLatencyPageCache.\n";
            } else {
                result_file << "LimitedMemLatencyPageCache achieves a " << (mem_hit_ratio - mira_hit_ratio) 
                           << "% higher hit ratio compared to MiraLatencyPageCache.\n";
            }
        }
        
        // Compare eviction patterns
        double eviction_ratio = (stats.mem_evictions > 0) ? 
            static_cast<double>(stats.mira_evictions) / stats.mem_evictions : 0;
            
        if (eviction_ratio > 10.0) {
            result_file << "Note: MiraLatencyPageCache performs significantly more evictions (" 
                       << stats.mira_evictions << " vs " << stats.mem_evictions 
                       << "), which may indicate more active cache management.\n";
        }
        
        result_file.close();
        
        // Print summary to console
        std::cout << "\n========== " << latency_name << " Latency Summary ==========\n";
        std::cout << std::fixed << std::setprecision(2);
        
        std::cout << "LimitedMemLatencyPageCache vs MiraLatencyPageCache performance comparison (median values):\n";
        std::cout << "  - Insert: " << mem_insert_median << " ms vs " << mira_insert_median 
                 << " ms (" << (insert_diff_pct > 0 ? "+" : "") << insert_diff_pct << "%)\n";
        
        std::cout << "  - Point Lookup: " << mem_lookup_median << " ms vs " << mira_lookup_median 
                 << " ms (" << (lookup_diff_pct > 0 ? "+" : "") << lookup_diff_pct << "%)\n";
        
        std::cout << "  - Range Scan: " << mem_range_median << " ms vs " << mira_range_median 
                 << " ms (" << (range_diff_pct > 0 ? "+" : "") << range_diff_pct << "%)\n";
        
        std::cout << "  - Mixed Workload: " << mem_mixed_median << " ms vs " << mira_mixed_median 
                 << " ms (" << (mixed_diff_pct > 0 ? "+" : "") << mixed_diff_pct << "%)\n";
        
        // Report cache hit ratios and evictions
        std::cout << "\nLimitedMemLatencyPageCache hit ratio: " << mem_hit_ratio << "%\n";
        std::cout << "LimitedMemLatencyPageCache evictions: " << stats.mem_evictions << "\n";
        std::cout << "MiraLatencyPageCache hit ratio: " << mira_hit_ratio << "%\n";
        std::cout << "MiraLatencyPageCache evictions: " << stats.mira_evictions << "\n";
        
        std::cout << "\nDetailed results written to: " << result_filename << "\n";
        std::cout << "=====================================\n";
    }
    
    // After all latency levels are benchmarked, generate CSV files for plotting
    if (generate_plots) {
        std::cout << "\nGenerating CSV files for plotting...\n";
        generate_csv_files(all_stats);
    }
}

int main(int argc, char* argv[]) {
    try {
        // Default values
        int num_iterations = 50;  // Default to 50 iterations for better statistics
        size_t page_size = 4096;  // Default page size (4KB)
        bool generate_plots = true;
        bool warm_up = true;
        
        // Parse command line arguments if provided
        if (argc > 1) {
            num_iterations = std::stoi(argv[1]);
        }
        
        if (argc > 2) {
            page_size = std::stoul(argv[2]);
        }
        
        if (argc > 3) {
            generate_plots = (std::stoi(argv[3]) != 0);
        }
        
        // Print banner
        std::cout << "==========================================\n";
        std::cout << "   Far Memory Latency Extended Benchmark   \n";
        std::cout << "==========================================\n";
        std::cout << "This benchmark will run " << num_iterations << " iterations\n";
        std::cout << "comparing MiraLatencyPageCache vs LimitedMemLatencyPageCache\n";
        std::cout << "with different simulated far memory latencies.\n\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Run benchmark with specified iterations
        run_multi_benchmark(num_iterations, page_size, generate_plots, warm_up);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        
        std::cout << "\nBenchmark completed in " << elapsed.count() << " seconds.\n";
        std::cout << "Results are available in the benchmark_results directory.\n";
        std::cout << "CSV files for plotting are available in the current directory.\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}