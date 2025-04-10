#ifndef _BPTREE_LIMITED_MEM_LATENCY_PAGE_CACHE_H_
#define _BPTREE_LIMITED_MEM_LATENCY_PAGE_CACHE_H_

#include "mira_page_cache_latency.h"
#include <list>
#include <mutex>

namespace bptree {

/**
 * A simplified page cache with size limits and latency simulation
 */
class LimitedMemLatencyPageCache : public AbstractPageCache {
public:
    LimitedMemLatencyPageCache(size_t page_size, LatencyLevel latency = LatencyLevel::NONE,
                             size_t max_cache_size_bytes = 4 * 1024 * 1024)
        : page_size(page_size), 
          latency(latency),
          max_pages(max_cache_size_bytes / page_size)
    {
        next_id.store(1);
        std::cout << "LimitedMemLatencyPageCache initialized with size limit: " 
                  << (max_cache_size_bytes / 1024) << " KB (" << max_pages << " pages)" << std::endl;
    }

    virtual Page* new_page(boost::upgrade_lock<Page>& lock) override
    {
        std::unique_lock<std::shared_mutex> guard(mutex);
        auto id = next_id.fetch_add(1);
        
        // Check if we need to evict pages
        while (page_map.size() >= max_pages && !lru_list.empty()) {
            // Evict the least recently used page
            PageID victim_id = lru_list.back();
            lru_list.pop_back();
            lru_map.erase(victim_id);
            page_map.erase(victim_id);
            evictions++;
        }
        
        // Create the new page
        auto* page = new Page(id, page_size);
        page_map[id] = std::unique_ptr<Page>(page);
        
        // Add to LRU tracking
        lru_list.push_front(id);
        lru_map[id] = lru_list.begin();
        
        lock = boost::upgrade_lock<Page>(*page);
        return page;
    }

    virtual Page* fetch_page(PageID id, boost::upgrade_lock<Page>& lock) override
    {
        // Add simulated latency for cache misses
        if (page_map.find(id) == page_map.end()) {
            add_simulated_latency();
            misses++;
            return nullptr;
        }
        
        std::unique_lock<std::shared_mutex> guard(mutex);
        auto it = page_map.find(id);
        if (it == page_map.end()) {
            misses++;
            return nullptr;
        }
        
        // Update LRU position
        auto lru_it = lru_map.find(id);
        if (lru_it != lru_map.end()) {
            lru_list.erase(lru_it->second);
            lru_list.push_front(id);
            lru_it->second = lru_list.begin();
        }
        
        hits++;
        lock = boost::upgrade_lock<Page>(*it->second);
        return it->second.get();
    }

    virtual void pin_page(Page* page, boost::upgrade_lock<Page>& lock) override {}
    virtual void unpin_page(Page* page, bool dirty, boost::upgrade_lock<Page>& lock) override {}
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
    size_t max_pages;
    LatencyLevel latency;
    std::atomic<PageID> next_id;
    std::shared_mutex mutex;
    std::unordered_map<PageID, std::unique_ptr<Page>> page_map;
    
    // LRU tracking
    std::list<PageID> lru_list;
    std::unordered_map<PageID, std::list<PageID>::iterator> lru_map;
    
    // Statistics
    std::atomic<uint64_t> hits{0};
    std::atomic<uint64_t> misses{0};
    std::atomic<uint64_t> evictions{0};
    
    // Add simulated latency for far memory access
    void add_simulated_latency() {
        switch (latency) {
            case LatencyLevel::LOW:
                std::this_thread::sleep_for(std::chrono::microseconds(100));
                break;
            case LatencyLevel::MEDIUM:
                std::this_thread::sleep_for(std::chrono::microseconds(500));
                break;
            case LatencyLevel::HIGH:
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                break;
            default:
                // No added latency
                break;
        }
    }
};

} // namespace bptree

#endif