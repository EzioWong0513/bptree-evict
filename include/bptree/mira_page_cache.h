#ifndef _BPTREE_MIRA_PAGE_CACHE_H_
#define _BPTREE_MIRA_PAGE_CACHE_H_

#include "page_cache.h"
#include "page.h"

#include <shared_mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <mutex>
#include <algorithm>
#include <string>
#include <functional>
#include <iostream>

namespace bptree {

// Forward declarations
class CacheSection;

/**
 * Enum representing different access patterns
 * as recognized in the Mira paper
 */
enum class AccessPattern {
    SEQUENTIAL,      // Sequential access pattern
    RANDOM,          // Random access with good locality
    STRIDE,          // Stride-based access
    UNKNOWN          // Unknown/default pattern
};

/**
 * CacheSection represents a portion of the cache dedicated to
 * a specific access pattern, mimicking Mira's section separation
 */
class CacheSection {
public:
    CacheSection(std::string name, AccessPattern pattern, size_t size, size_t line_size)
        : name(name), pattern(pattern), size(size), line_size(line_size) {}

    // Get cache section info
    std::string get_name() const { return name; }
    AccessPattern get_pattern() const { return pattern; }
    size_t get_size() const { return size; }
    size_t get_line_size() const { return line_size; }

    // Add a page to this section
    void add_page(PageID id) {
        pages.push_front(id);
        page_it_map[id] = pages.begin();
    }

    // Remove a page from this section
    void remove_page(PageID id) {
        auto it = page_it_map.find(id);
        if (it != page_it_map.end()) {
            pages.erase(it->second);
            page_it_map.erase(it);
        }
    }

    // Mark a page as evictable (Mira's program-guided eviction hint)
    void mark_evictable(PageID id) {
        evictable_pages.insert(id);
    }

    // Check if a page is marked as evictable
    bool is_evictable(PageID id) const {
        return evictable_pages.find(id) != evictable_pages.end();
    }

    // Select victim page for eviction based on section's access pattern
    bool select_victim(PageID& victim_id) {
        if (pages.empty())
            return false;

        // First try to find evictable pages (Mira's hint-based eviction)
        for (auto it = pages.rbegin(); it != pages.rend(); ++it) {
            if (is_evictable(*it)) {
                victim_id = *it;
                pages.erase(std::next(it).base());
                page_it_map.erase(victim_id);
                evictable_pages.erase(victim_id);
                return true;
            }
        }

        // If no evictable pages, use pattern-specific strategy
        switch (pattern) {
            case AccessPattern::SEQUENTIAL:
                // For sequential, oldest pages are least likely to be reused
                victim_id = pages.back();
                pages.pop_back();
                page_it_map.erase(victim_id);
                return true;
            
            case AccessPattern::RANDOM:
                // For random with locality, use LRU
                victim_id = pages.back();
                pages.pop_back();
                page_it_map.erase(victim_id);
                return true;
                
            case AccessPattern::STRIDE:
                // For stride patterns, consider stride distance
                // (simplified to LRU in this implementation)
                victim_id = pages.back();
                pages.pop_back();
                page_it_map.erase(victim_id);
                return true;
                
            default:
                // Default to LRU
                victim_id = pages.back();
                pages.pop_back();
                page_it_map.erase(victim_id);
                return true;
        }
    }

    // Access a page, updating its position in the appropriate data structure
    void touch_page(PageID id) {
        auto it = page_it_map.find(id);
        if (it != page_it_map.end()) {
            pages.erase(it->second);
            pages.push_front(id);
            it->second = pages.begin();
        }
    }

    // Get current page count
    size_t page_count() const {
        return pages.size();
    }

private:
    std::string name;
    AccessPattern pattern;
    size_t size;         // Size in bytes
    size_t line_size;    // Cache line size in bytes
    
    // List of pages in this section (front is most recently used)
    std::list<PageID> pages;
    
    // Map for O(1) lookup of list iterators
    std::unordered_map<PageID, std::list<PageID>::iterator> page_it_map;
    
    // Set of pages marked as evictable by program hints
    std::unordered_set<PageID> evictable_pages;
};

/**
 * MiraPageCache implements a Mira-inspired page cache with
 * multiple sections for different access patterns and
 * program-guided eviction hints
 */
class MiraPageCache : public AbstractPageCache {
public:
    MiraPageCache(size_t page_size, size_t total_cache_size = 4 * 1024 * 1024)
        : page_size(page_size), total_cache_size(total_cache_size)
    {
        next_id.store(1);
        
        // Initialize cache sections with different configurations
        // mimicking Mira's section separation
        size_t seq_size = total_cache_size * 0.2;  // 20% for sequential access
        size_t random_size = total_cache_size * 0.6;  // 60% for random access
        size_t stride_size = total_cache_size * 0.2;  // 20% for stride access
        
        // Create sections with appropriate patterns and line sizes
        sections.emplace_back("sequential", AccessPattern::SEQUENTIAL, seq_size, page_size * 4);
        sections.emplace_back("random", AccessPattern::RANDOM, random_size, page_size);
        sections.emplace_back("stride", AccessPattern::STRIDE, stride_size, page_size * 2);
        
        // Default section for pages with unknown pattern
        default_section = std::make_unique<CacheSection>("default", AccessPattern::UNKNOWN, 
                                                        total_cache_size, page_size);
                                                        
        std::cout << "MiraPageCache initialized with " << sections.size() 
                  << " sections and a default section" << std::endl;
    }

    virtual Page* new_page(boost::upgrade_lock<Page>& lock) override
    {
        auto id = get_next_id();
        std::unique_lock<std::shared_mutex> guard(mutex);
        page_map[id] = std::make_unique<Page>(id, page_size);
        Page* page = page_map[id].get();
        
        // Assign to default section initially
        page_section_map[id] = -1; // -1 means default section
        default_section->add_page(id);
        
        lock = boost::upgrade_lock<Page>(*page);
        return page;
    }

    virtual Page* fetch_page(PageID id, boost::upgrade_lock<Page>& lock) override
    {
        std::shared_lock<std::shared_mutex> guard(mutex);
        auto it = page_map.find(id);
        if (it == page_map.end()) {
            // Page not in cache
            misses++;
            return nullptr;
        }
        
        // Page found in cache
        hits++;
        
        lock = boost::upgrade_lock<Page>(*it->second);
        
        // Update LRU position in its section
        update_lru(id);
        
        return it->second.get();
    }

    virtual void pin_page(Page* page, boost::upgrade_lock<Page>& lock) override 
    {
        // When a page is pinned, remove it from the eviction candidates
        PageID id = page->get_id();
        std::lock_guard<std::mutex> guard(section_mutex);
        
        int section_idx = page_section_map[id];
        if (section_idx >= 0) {
            sections[section_idx].touch_page(id);
        } else if (section_idx == -1) {
            default_section->touch_page(id);
        }
    }
    
    virtual void unpin_page(Page* page, bool dirty, boost::upgrade_lock<Page>& lock) override 
    {
        // Add the page back to eviction candidates when unpinned
        // (No-op for our simplified implementation)
    }

    // Implement Mira's program-guided eviction hint
    void mark_evictable(PageID id) 
    {
        std::lock_guard<std::mutex> guard(section_mutex);
        
        int section_idx = page_section_map[id];
        if (section_idx >= 0) {
            sections[section_idx].mark_evictable(id);
        } else if (section_idx == -1) {
            default_section->mark_evictable(id);
        }
    }

    // Assign a page to a specific cache section based on its access pattern
    void assign_section(PageID id, AccessPattern pattern) 
    {
        std::lock_guard<std::mutex> guard(section_mutex);
        
        // First remove from current section
        int current_section = page_section_map[id];
        if (current_section >= 0) {
            sections[current_section].remove_page(id);
        } else if (current_section == -1) {
            default_section->remove_page(id);
        }
        
        // Find appropriate section for the pattern
        for (size_t i = 0; i < sections.size(); i++) {
            if (sections[i].get_pattern() == pattern) {
                sections[i].add_page(id);
                page_section_map[id] = i;
                return;
            }
        }
        
        // If no matching section, use default
        default_section->add_page(id);
        page_section_map[id] = -1;
    }

    virtual void flush_page(Page* page, boost::upgrade_lock<Page>& lock) override {}
    virtual void flush_all_pages() override {}

    virtual size_t size() const override { return page_map.size(); }
    virtual size_t get_page_size() const override { return page_size; }
    
    // Stats methods for benchmarking
    void reset_stats() {
        hits = 0;
        misses = 0;
        evictions = 0;
    }
    
    uint64_t get_hits() const { return hits; }
    uint64_t get_misses() const { return misses; }
    uint64_t get_evictions() const { return evictions; }

private:
    size_t page_size;
    size_t total_cache_size;
    std::atomic<PageID> next_id;
    std::shared_mutex mutex;
    std::mutex section_mutex;
    
    // Page storage
    std::unordered_map<PageID, std::unique_ptr<Page>> page_map;
    
    // Cache sections for different access patterns
    std::vector<CacheSection> sections;
    std::unique_ptr<CacheSection> default_section;
    
    // Maps page IDs to section indexes (-1 for default section)
    std::unordered_map<PageID, int> page_section_map;
    
    // Statistics for benchmarking
    std::atomic<uint64_t> hits{0};
    std::atomic<uint64_t> misses{0};
    std::atomic<uint64_t> evictions{0};

    PageID get_next_id() { return next_id++; }

    // Update a page's position in its section's LRU list
    void update_lru(PageID id) {
        std::lock_guard<std::mutex> guard(section_mutex);
        
        int section_idx = page_section_map[id];
        if (section_idx >= 0) {
            sections[section_idx].touch_page(id);
        } else if (section_idx == -1) {
            default_section->touch_page(id);
        }
    }
    
    // Find a victim page for eviction
    bool find_victim(PageID& victim_id) {
        // First try sections with evictable pages
        for (auto& section : sections) {
            PageID id;
            if (section.select_victim(id)) {
                victim_id = id;
                evictions++;
                return true;
            }
        }
        
        // Then try default section
        if (default_section->select_victim(victim_id)) {
            evictions++;
            return true;
        }
        
        // If no victims found, return false
        return false;
    }
};

} // namespace bptree

#endif