#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Set up the style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set(font_scale=1.3)
colors = sns.color_palette("muted")

# Create results directory if it doesn't exist
os.makedirs("benchmark_plots", exist_ok=True)

# Function to generate plots for an operation
def generate_operation_plots(operation_name, file_name):
    """Generate performance comparison plots for a specific operation"""
    # Load the data
    try:
        df = pd.read_csv(f"{file_name}.csv")
    except FileNotFoundError:
        print(f"Warning: {file_name}.csv not found, skipping {operation_name} plots")
        return
    
    print(f"Generating plots for {operation_name}...")
    
    # Plot 1: Performance comparison across latency levels (boxplot)
    plt.figure(figsize=(12, 8))
    sns.boxplot(x="LatencyLevel", y="Baseline", data=df, color=colors[0])
    sns.boxplot(x="LatencyLevel", y="MiraCache", data=df, color=colors[1])
    plt.title(f"{operation_name} Performance Comparison Across Latency Levels")
    plt.xlabel("Latency Level")
    plt.ylabel("Time (ms)")
    plt.legend(["BaselineCache", "MiraLatencyPageCache"])
    plt.savefig(f"benchmark_plots/{operation_name.lower().replace(' ', '_')}_boxplot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Performance difference percentage across latency levels (violin plot)
    plt.figure(figsize=(12, 8))
    sns.violinplot(x="LatencyLevel", y="Difference(%)", data=df, palette="muted")
    plt.title(f"{operation_name} Performance Difference % (MiraCache vs Baseline)")
    plt.xlabel("Latency Level")
    plt.ylabel("Difference %")
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)  # Add horizontal line at y=0
    plt.savefig(f"benchmark_plots/{operation_name.lower().replace(' ', '_')}_difference.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Performance comparison over iterations for each latency level
    latency_levels = df["LatencyLevel"].unique()
    
    for latency in latency_levels:
        latency_df = df[df["LatencyLevel"] == latency]
        
        plt.figure(figsize=(12, 6))
        plt.plot(latency_df["Iteration"], latency_df["Baseline"], label="BaselineCache", marker='o')
        plt.plot(latency_df["Iteration"], latency_df["MiraCache"], label="MiraLatencyPageCache", marker='x')
        plt.title(f"{operation_name} Performance Over Iterations ({latency} latency)")
        plt.xlabel("Iteration")
        plt.ylabel("Time (ms)")
        plt.legend()
        plt.grid(True)
        plt.savefig(f"benchmark_plots/{operation_name.lower().replace(' ', '_')}_{latency}_iterations.png", dpi=300, bbox_inches='tight')
        plt.close()

# Function to generate hit ratio plots
def generate_hit_ratio_plots():
    """Generate plots for cache hit ratios"""
    try:
        df = pd.read_csv("hit_ratios.csv")
    except FileNotFoundError:
        print("Warning: hit_ratios.csv not found, skipping hit ratio plots")
        return
    
    print("Generating hit ratio plots...")
    
    # Plot hit ratios by operation and latency level
    plt.figure(figsize=(15, 8))
    
    # Reshape data for easier plotting
    melted_df = pd.melt(df, id_vars=["Operation", "LatencyLevel"],
                       value_vars=["BaselineHitRatio", "MiraCacheHitRatio"],
                       var_name="CacheType", value_name="HitRatio")
    
    # Create grouped bar plot
    ax = sns.barplot(x="Operation", y="HitRatio", hue="CacheType", data=melted_df, 
                    palette=[colors[0], colors[1]], ci=None)
    
    # Add labels and title
    plt.title("Cache Hit Ratio Comparison")
    plt.xlabel("Operation")
    plt.ylabel("Hit Ratio (%)")
    plt.ylim(0, 100)  # Set y-axis limits
    
    # Update legend labels
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, ["BaselineCache", "MiraLatencyPageCache"])
    
    # Add value labels on top of bars
    for p in ax.patches:
        ax.annotate(f'{p.get_height():.1f}%', 
                   (p.get_x() + p.get_width() / 2., p.get_height()), 
                   ha = 'center', va = 'bottom',
                   xytext = (0, 5), textcoords = 'offset points')
    
    # Save the plot
    plt.savefig("benchmark_plots/hit_ratio_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot hit ratio by latency level for each operation
    operations = df["Operation"].unique()
    
    for operation in operations:
        op_df = df[df["Operation"] == operation]
        
        plt.figure(figsize=(10, 6))
        
        # Bar plot for each latency level
        x = np.arange(len(op_df))
        width = 0.35
        
        plt.bar(x - width/2, op_df["BaselineHitRatio"], width, label="BaselineCache", color=colors[0])
        plt.bar(x + width/2, op_df["MiraCacheHitRatio"], width, label="MiraLatencyPageCache", color=colors[1])
        
        plt.title(f"Hit Ratio Comparison for {operation} Operation")
        plt.xlabel("Latency Level")
        plt.ylabel("Hit Ratio (%)")
        plt.xticks(x, op_df["LatencyLevel"])
        plt.ylim(0, 100)
        plt.legend()
        
        # Add value labels on bars
        for i, v in enumerate(op_df["BaselineHitRatio"]):
            plt.text(i - width/2, v + 1, f"{v:.1f}%", ha='center')
        
        for i, v in enumerate(op_df["MiraCacheHitRatio"]):
            plt.text(i + width/2, v + 1, f"{v:.1f}%", ha='center')
        
        plt.savefig(f"benchmark_plots/hit_ratio_{operation.lower()}.png", dpi=300, bbox_inches='tight')
        plt.close()

# Function to generate eviction comparison plots
def generate_eviction_plots():
    """Generate plots for cache eviction patterns"""
    try:
        df = pd.read_csv("eviction_stats.csv")
    except FileNotFoundError:
        print("Warning: eviction_stats.csv not found, skipping eviction plots")
        return
    
    print("Generating eviction pattern plots...")
    
    # Plot 1: Eviction count comparison
    plt.figure(figsize=(12, 7))
    
    x = np.arange(len(df))
    width = 0.35
    
    plt.bar(x - width/2, df["BaselineEvictions"], width, label="BaselineCache", color=colors[0])
    plt.bar(x + width/2, df["MiraCacheEvictions"], width, label="MiraLatencyPageCache", color=colors[1])
    
    plt.title("Cache Eviction Count Comparison")
    plt.xlabel("Latency Level")
    plt.ylabel("Number of Evictions")
    plt.xticks(x, df["LatencyLevel"])
    plt.legend()
    
    # Add value labels on bars
    for i, v in enumerate(df["BaselineEvictions"]):
        plt.text(i - width/2, v + 0.5, f"{v}", ha='center')
    
    for i, v in enumerate(df["MiraCacheEvictions"]):
        plt.text(i + width/2, v + 0.5, f"{v}", ha='center')
    
    plt.savefig("benchmark_plots/eviction_count_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Eviction ratio (MiraCache/Baseline)
    plt.figure(figsize=(10, 6))
    plt.bar(df["LatencyLevel"], df["EvictionRatio"], color=colors[2])
    plt.title("Eviction Ratio (MiraCache/Baseline)")
    plt.xlabel("Latency Level")
    plt.ylabel("Ratio")
    
    # Add value labels on bars
    for i, v in enumerate(df["EvictionRatio"]):
        plt.text(i, v + 0.1, f"{v:.2f}x", ha='center')
    
    plt.savefig("benchmark_plots/eviction_ratio.png", dpi=300, bbox_inches='tight')
    plt.close()

# Function to generate summary comparison plots
def generate_summary_plots():
    """Generate summary comparison plots"""
    try:
        df = pd.read_csv("summary_stats.csv")
    except FileNotFoundError:
        print("Warning: summary_stats.csv not found, skipping summary plots")
        return
    
    print("Generating summary comparison plots...")
    
    # Plot 1: Performance improvement by operation and latency level
    plt.figure(figsize=(15, 10))
    
    # Filter for necessary columns
    plot_df = df[["Operation", "LatencyLevel", "DifferenceMedian(%)"]]
    
    # Create the bar plot
    ax = sns.barplot(x="Operation", y="DifferenceMedian(%)", hue="LatencyLevel", data=plot_df)
    
    plt.title("Performance Difference by Operation and Latency Level")
    plt.xlabel("Operation")
    plt.ylabel("Difference % (MiraCache vs Baseline)")
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)  # Add horizontal line at y=0
    
    # Add value labels on bars
    for p in ax.patches:
        value = f"{p.get_height():.1f}%"
        ax.annotate(value, 
                   (p.get_x() + p.get_width() / 2., p.get_height()),
                   ha = 'center', va = 'bottom' if p.get_height() > 0 else 'top',
                   xytext = (0, 5 if p.get_height() > 0 else -10), 
                   textcoords = 'offset points')
    
    plt.savefig("benchmark_plots/performance_difference_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Performance comparison across operations for each latency level
    latency_levels = df["LatencyLevel"].unique()
    
    for latency in latency_levels:
        latency_df = df[df["LatencyLevel"] == latency]
        
        plt.figure(figsize=(12, 8))
        
        # Reshape data for easier plotting
        plot_data = pd.melt(latency_df, 
                           id_vars=["Operation"],
                           value_vars=["BaselineMedian", "MiraCacheMedian"],
                           var_name="CacheType", value_name="Time (ms)")
        
        # Create the bar plot
        ax = sns.barplot(x="Operation", y="Time (ms)", hue="CacheType", data=plot_data)
        
        plt.title(f"Performance Comparison by Operation ({latency} latency)")
        plt.xlabel("Operation")
        plt.ylabel("Time (ms)")
        
        # Update legend labels
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, ["BaselineCache", "MiraLatencyPageCache"])
        
        # Add value labels on bars
        for p in ax.patches:
            ax.annotate(f'{p.get_height():.1f}', 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'bottom',
                       xytext = (0, 5), textcoords = 'offset points')
        
        plt.savefig(f"benchmark_plots/performance_comparison_{latency}.png", dpi=300, bbox_inches='tight')
        plt.close()

# Function to generate combined performance plots
def generate_combined_performance_plot():
    """Generate a combined performance comparison plot"""
    try:
        summary_df = pd.read_csv("summary_stats.csv")
    except FileNotFoundError:
        print("Warning: summary_stats.csv not found, skipping combined performance plot")
        return
    
    print("Generating combined performance plot...")
    
    # Find the most important operation - Point Lookup
    lookup_df = summary_df[summary_df["Operation"] == "Lookup"]
    
    # Reshape data for easier plotting
    plot_data = pd.melt(lookup_df, 
                       id_vars=["LatencyLevel"],
                       value_vars=["BaselineMedian", "MiraCacheMedian"],
                       var_name="CacheType", value_name="Time (ms)")
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x="LatencyLevel", y="Time (ms)", hue="CacheType", data=plot_data)
    
    plt.title("Point Lookup Performance Comparison")
    plt.xlabel("Latency Level")
    plt.ylabel("Time (ms)")
    
    # Update legend labels
    plt.legend(["BaselineCache", "MiraLatencyPageCache"])
    
    plt.savefig("benchmark_plots/combined_lookup_performance.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a difference percentage plot across latency levels
    plt.figure(figsize=(10, 6))
    sns.barplot(x="LatencyLevel", y="DifferenceMedian(%)", data=lookup_df, color=colors[2])
    
    plt.title("Point Lookup Performance Difference % (MiraCache vs Baseline)")
    plt.xlabel("Latency Level")
    plt.ylabel("Difference %")
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.3)  # Add horizontal line at y=0
    
    # Add value labels on bars
    for i, v in enumerate(lookup_df["DifferenceMedian(%)"]):
        color = "green" if v < 0 else "red"
        plt.text(i, v + 0.1 if v >= 0 else v - 0.5, f"{v:.2f}%", ha='center', color=color)
    
    plt.savefig("benchmark_plots/combined_lookup_difference.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    print("Generating benchmark visualization plots...")
    
    # Generate operation-specific plots
    generate_operation_plots("Point Lookup", "lookup_times")
    generate_operation_plots("Range Scan", "range_scan_times")
    generate_operation_plots("Mixed Workload", "mixed_workload_times")
    
    # Generate hit ratio plots
    generate_hit_ratio_plots()
    
    # Generate eviction plots
    generate_eviction_plots()
    
    # Generate summary plots
    generate_summary_plots()
    
    # Generate combined performance plot
    generate_combined_performance_plot()
    
    print("All plots generated and saved in the 'benchmark_plots' directory.")

if __name__ == "__main__":
    main()