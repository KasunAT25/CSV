#ifndef __LIPP_H__
#define __LIPP_H__

#include "lipp_base.h"
#include <stdint.h>
#include <math.h>
#include <limits>
#include <cstdio>
#include <stack>
#include <vector>
#include <cstring>
#include <sstream>

typedef uint8_t bitmap_t;
#define BITMAP_WIDTH (sizeof(bitmap_t) * 8)
#define BITMAP_SIZE(num_items) (((num_items) + BITMAP_WIDTH - 1) / BITMAP_WIDTH)
#define BITMAP_GET(bitmap, pos) (((bitmap)[(pos) / BITMAP_WIDTH] >> ((pos) % BITMAP_WIDTH)) & 1)
#define BITMAP_SET(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] |= 1 << ((pos) % BITMAP_WIDTH))
#define BITMAP_CLEAR(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] &= ~bitmap_t(1 << ((pos) % BITMAP_WIDTH)))
#define BITMAP_NEXT_1(bitmap_item) __builtin_ctz((bitmap_item))
// TO REMOVE POI_BITMAP
// Delete all instances of poi_bitmap amd change else if to else 
// runtime assert
//I changed the struct Node to public and put it at the top as well as Items
#define RT_ASSERT(expr) \
{ \
    if (!(expr)) { \
        fprintf(stderr, "RT_ASSERT Error at %s:%d, `%s`\n", __FILE__, __LINE__, #expr); \
        exit(0); \
    } \
}

#define COLLECT_TIME 0

#if COLLECT_TIME
#include <chrono>
#endif

template<class T, class P, bool USE_FMCD = true>
class LIPP
{
public:
        static_assert(std::is_arithmetic<T>::value, "LIPP key type must be numeric.");

        inline int compute_gap_count(int size) {
            if (size >= 1000000) return 1;
            if (size >= 100000) return 2;
            return 5;
        }
        //struct Node;
    struct Node;
    struct Item
        {
            union {
                struct {
                    T key;
                    P value;
                } data;
                Node* child;
            } comp;
        };
    struct Node
        {
            int is_two; // is special node for only two keys
            int build_size; // tree size (include sub nodes) when node created
            int size; // current tree size (include sub nodes)
            int fixed; // fixed node will not trigger rebuild
            int num_inserts, num_insert_to_data;
            int num_items; // size of items
            LinearModel<T> model;
            Item* items;
            bitmap_t* none_bitmap; // 1 means None, 0 means Data or Child
            bitmap_t* child_bitmap; // 1 means Child. will always be 0 when none_bitmap is 1

            bitmap_t* poi_bitmap; // 1 means poisoned. 0 not
            int poi_size = 0;

            //Need a new one

            int level = 1; //I added this to make the level easier
        };

        inline int PREDICT_POS(Node* node, T key) const {
            double v = node->model.predict_double(key);
            if (v > std::numeric_limits<int>::max() / 2) {
                return node->num_items - 1;
            }
            if (v < 0) {
                return 0;
            }
            return std::min(node->num_items - 1, static_cast<int>(v));
        }

        static void remove_last_bit(bitmap_t& bitmap_item) {
            bitmap_item -= 1 << BITMAP_NEXT_1(bitmap_item);
        }

        const double BUILD_LR_REMAIN;
        const bool QUIET;

        struct {
            long long fmcd_success_times = 0;
            long long fmcd_broken_times = 0;
            #if COLLECT_TIME
            double time_scan_and_destory_tree = 0;
            double time_build_tree_bulk = 0;
            #endif
        } stats;

public:
    typedef std::pair<T, P> V;

    LIPP(double BUILD_LR_REMAIN = 0, bool QUIET = true)
        : BUILD_LR_REMAIN(BUILD_LR_REMAIN), QUIET(QUIET) {
        {
            std::vector<Node*> nodes;
            for (int _ = 0; _ < 1e7; _ ++) {
                Node* node = build_tree_two(T(0), P(), T(1), P());
                nodes.push_back(node);
            }
            for (auto node : nodes) {
                destroy_tree(node);
            }
            if (!QUIET) {
                printf("initial memory pool size = %lu\n", pending_two.size());
            }
        }
        if (USE_FMCD && !QUIET) {
            printf("enable FMCD\n");
        }

        root = build_tree_none();
    }
    ~LIPP() {
        destroy_tree(root);
        root = NULL;
        destory_pending();
    }

    void insert(const V& v) {
        insert(v.first, v.second);
    }
    void insert(const T& key, const P& value) {
        //I changed from insert_tree to insert_tree2
        root = insert_tree2(root, key, value);
    }
    P at(const T& key, bool skip_existence_check = true) const {
        Node* node = root;

        while (true) {
            int pos = PREDICT_POS(node, key);
            if (BITMAP_GET(node->child_bitmap, pos) == 1) {
                node = node->items[pos].comp.child;
            } else {
                if (skip_existence_check) {
                    return node->items[pos].comp.data.value;
                } else {
                    if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                        RT_ASSERT(false);
                    } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
                        RT_ASSERT(node->items[pos].comp.data.key == key);
                        return node->items[pos].comp.data.value;
                    }
                }
            }
        }
    }
    bool exists(const T& key) const {
        Node* node = root;
        while (true) {
            int pos = PREDICT_POS(node, key);
            if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                return false;
            } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
                return node->items[pos].comp.data.key == key;
            } else {
                node = node->items[pos].comp.child;
            }
        }
    }
    void bulk_load(const V* vs, int num_keys) {
        if (num_keys == 0) {
            destroy_tree(root);
            root = build_tree_none();
            return;
        }
        if (num_keys == 1) {
            destroy_tree(root);
            root = build_tree_none();
            insert(vs[0]);
            return;
        }
        if (num_keys == 2) {
            destroy_tree(root);
            root = build_tree_two(vs[0].first, vs[0].second, vs[1].first, vs[1].second);
            return;
        }

        RT_ASSERT(num_keys > 2);
        for (int i = 1; i < num_keys; i ++) {
            RT_ASSERT(vs[i].first > vs[i-1].first);
        }

        T* keys = new T[num_keys];
        P* values = new P[num_keys];
        for (int i = 0; i < num_keys; i ++) {
            keys[i] = vs[i].first;
            values[i] = vs[i].second;
        }
        destroy_tree(root);
        root = build_tree_bulk(keys, values, num_keys);
        delete[] keys;
        delete[] values;
    }

    void show() const {
        printf("============= SHOW LIPP ================\n");

        std::stack<Node*> s;
        s.push(root);
        while (!s.empty()) {
            Node* node = s.top(); s.pop();

            printf("Node(%p, a = %lf, b = %lf, num_items = %d)", node, node->model.a, node->model.b, node->num_items);
            printf("[");
            int first = 1;
            for (int i = 0; i < node->num_items; i ++) {
                if (!first) {
                    printf(", ");
                }
                first = 0;
                if (BITMAP_GET(node->none_bitmap, i) == 1) {
                    printf("None");
                } else if (BITMAP_GET(node->child_bitmap, i) == 0) {
                    std::stringstream s;
                    s << node->items[i].comp.data.key;
                    printf("Key(%s)", s.str().c_str());
                } else {
                    printf("Child(%p)", node->items[i].comp.child);
                    s.push(node->items[i].comp.child);
                }
            }
            printf("]\n");
        }
    }
    void print_depth() const {
        std::stack<Node*> s;
        std::stack<int> d;
        s.push(root);
        d.push(1);

        int max_depth = 1;
        int sum_depth = 0, sum_nodes = 0;
        while (!s.empty()) {
            Node* node = s.top(); s.pop();
            int depth = d.top(); d.pop();
            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->child_bitmap, i) == 1) {
                    s.push(node->items[i].comp.child);
                    d.push(depth + 1);
                } else if (BITMAP_GET(node->none_bitmap, i) != 1) {
                    max_depth = std::max(max_depth, depth);
                    sum_depth += depth;
                    sum_nodes ++;
                }
            }
        }

        printf("max_depth = %d, avg_depth = %.2lf\n", max_depth, double(sum_depth) / double(sum_nodes));
    }
    //To get all nodes in the index and set their level.
    void scan_nodes(std::vector<Node*> *nodes) {
        //Starts from the root
        std::stack<Node*> s;
        std::stack<int> d;

        //Push the root
        s.push(root);
        //this is the depth
        d.push(1);

        int max_depth = 1;
        int sum_depth = 0, sum_nodes = 0;
        while (!s.empty()) {
            Node* node = s.top(); s.pop();
            int depth = d.top(); d.pop();
            node->level = depth;
            nodes->push_back(node);

            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->child_bitmap, i) == 1) {
                    s.push(node->items[i].comp.child);
                    d.push(depth + 1);
                } else if (BITMAP_GET(node->none_bitmap, i) != 1) {
                    max_depth = std::max(max_depth, depth);
                    sum_depth += depth;
                    sum_nodes ++;
                }
            }
        }

        printf("max_depth = %d, avg_depth = %.2lf\n", max_depth, double(sum_depth) / double(sum_nodes));
    }
    void verify() const {
        std::stack<Node*> s;
        s.push(root);

        while (!s.empty()) {
            Node* node = s.top(); s.pop();
            int sum_size = 0;
            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->child_bitmap, i) == 1) {
                    s.push(node->items[i].comp.child);
                    sum_size += node->items[i].comp.child->size;
                } else if (BITMAP_GET(node->none_bitmap, i) != 1) {
                    sum_size ++;
                }
            }
            RT_ASSERT(sum_size == node->size);
        }
    }
    void print_stats() const {
        printf("======== Stats ===========\n");
        if (USE_FMCD) {
            printf("\t fmcd_success_times = %lld\n", stats.fmcd_success_times);
            printf("\t fmcd_broken_times = %lld\n", stats.fmcd_broken_times);
        }
        #if COLLECT_TIME
        printf("\t time_scan_and_destory_tree = %lf\n", stats.time_scan_and_destory_tree);
        printf("\t time_build_tree_bulk = %lf\n", stats.time_build_tree_bulk);
        #endif
    }
    size_t index_size(bool total=false, bool ignore_child=true) const {
        std::stack<Node*> s;
        s.push(root);
    
        size_t size = 0;
        while (!s.empty()) {
            Node* node = s.top(); s.pop();
            bool has_child = false;
            if(ignore_child == false) {
                size += sizeof(*node);
            }
            for (int i = 0; i < node->num_items; i ++) {
                if (ignore_child == true) {
                    size += sizeof(Item);
                    has_child = true;
                } else {
                    if (total) size += sizeof(Item);
                }
                if (BITMAP_GET(node->child_bitmap, i) == 1) {
                    if (!total) size += sizeof(Item);
                    s.push(node->items[i].comp.child);
                }
            }
            if (ignore_child == true && has_child) {
                size += sizeof(*node);
            }
        }
        return size;
    }

public:
    //struct Node;
    // struct Item
    // {
    //     union {
    //         struct {
    //             T key;
    //             P value;
    //         } data;
    //         Node* child;
    //     } comp;
    // };
    // struct Node
    // {
    //     int is_two; // is special node for only two keys
    //     int build_size; // tree size (include sub nodes) when node created
    //     int size; // current tree size (include sub nodes)
    //     int fixed; // fixed node will not trigger rebuild
    //     int num_inserts, num_insert_to_data;
    //     int num_items; // size of items
    //     LinearModel<T> model;
    //     Item* items;
    //     bitmap_t* none_bitmap; // 1 means None, 0 means Data or Child
    //     bitmap_t* child_bitmap; // 1 means Child. will always be 0 when none_bitmap is 1
    // };

    Node* root;
    std::stack<Node*> pending_two;

    std::allocator<Node> node_allocator;
    Node* new_nodes(int n)
    {
        Node* p = node_allocator.allocate(n);
        RT_ASSERT(p != NULL && p != (Node*)(-1));
        return p;
    }
    void delete_nodes(Node* p, int n)
    {
        node_allocator.deallocate(p, n);
    }

    std::allocator<Item> item_allocator;
    Item* new_items(int n)
    {
        Item* p = item_allocator.allocate(n);
        RT_ASSERT(p != NULL && p != (Item*)(-1));
        return p;
    }
    void delete_items(Item* p, int n)
    {
        item_allocator.deallocate(p, n);
    }

    std::allocator<bitmap_t> bitmap_allocator;
    bitmap_t* new_bitmap(int n)
    {
        bitmap_t* p = bitmap_allocator.allocate(n);
        RT_ASSERT(p != NULL && p != (bitmap_t*)(-1));
        return p;
    }
    void delete_bitmap(bitmap_t* p, int n)
    {
        bitmap_allocator.deallocate(p, n);
    }

    /// build an empty tree
    Node* build_tree_none()
    {
        Node* node = new_nodes(1);
        node->is_two = 0;
        node->build_size = 0;
        node->size = 0;
        node->fixed = 0;
        node->num_inserts = node->num_insert_to_data = 0;
        node->num_items = 1;
        node->model.a = node->model.b = 0;
        node->items = new_items(1);
        node->none_bitmap = new_bitmap(1);
        node->none_bitmap[0] = 0;
        BITMAP_SET(node->none_bitmap, 0);
        node->child_bitmap = new_bitmap(1);
        node->child_bitmap[0] = 0;

        node->poi_bitmap = new_bitmap(1);
        node->poi_bitmap[0] = 0;

        return node;
    }
    /// build a tree with two keys
    Node* build_tree_two(T key1, P value1, T key2, P value2)
    {
        if (key1 > key2) {
            std::swap(key1, key2);
            std::swap(value1, value2);
        }
        // std::cout << key1 << std::endl;
        // std::cout << key2 << std::endl;

        RT_ASSERT(key1 < key2);
        static_assert(BITMAP_WIDTH == 8);

        Node* node = NULL;
        if (pending_two.empty()) {
            node = new_nodes(1);
            node->is_two = 1;
            node->build_size = 2;
            node->size = 2;
            node->fixed = 0;
            node->num_inserts = node->num_insert_to_data = 0;

            node->num_items = 8;
            node->items = new_items(node->num_items);
            node->none_bitmap = new_bitmap(1);
            node->child_bitmap = new_bitmap(1);
            node->none_bitmap[0] = 0xff;
            node->child_bitmap[0] = 0;

            node->poi_bitmap = new_bitmap(1);
            node->poi_bitmap[0] = 0;
        } else {
            node = pending_two.top(); pending_two.pop();
        }

        const long double mid1_key = key1;
        const long double mid2_key = key2;

        const double mid1_target = node->num_items / 3;
        const double mid2_target = node->num_items * 2 / 3;

        node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
        node->model.b = mid1_target - node->model.a * mid1_key;
        RT_ASSERT(isfinite(node->model.a));
        RT_ASSERT(isfinite(node->model.b));

        { // insert key1&value1
            int pos = PREDICT_POS(node, key1);
            RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
            BITMAP_CLEAR(node->none_bitmap, pos);
            node->items[pos].comp.data.key = key1;
            node->items[pos].comp.data.value = value1;
        }
        { // insert key2&value2
            int pos = PREDICT_POS(node, key2);
            RT_ASSERT(BITMAP_GET(node->none_bitmap, pos) == 1);
            BITMAP_CLEAR(node->none_bitmap, pos);
            node->items[pos].comp.data.key = key2;
            node->items[pos].comp.data.value = value2;
        }

        return node;
    }
    /// bulk build, _keys must be sorted in asc order.
    Node* build_tree_bulk(T* _keys, P* _values, int _size)
    {
        if (USE_FMCD) {
            return build_tree_bulk_fmcd(_keys, _values, _size);
        } else {
            return build_tree_bulk_fast(_keys, _values, _size);
        }
    }

    ///For poisoning without adding the poisoning. bulk build, _keys must be sorted in asc order.
    Node* poi_build_tree_bulk(T* _keys, P* _values, int _size, T* _keys_real, int _size_real)
    {
        if (USE_FMCD) {
            return poi_build_tree_bulk_fmcd(_keys, _values, _size,  _keys_real, _size_real);
        } else {
            return build_tree_bulk_fast(_keys, _values, _size);
        }
    }
    /// bulk build, _keys must be sorted in asc order.
    /// split keys into three parts at each node.
    Node* build_tree_bulk_fast(T* _keys, P* _values, int _size)
    {
        RT_ASSERT(_size > 1);

        typedef struct {
            int begin;
            int end;
            int level; // top level = 1
            Node* node;
        } Segment;
        std::stack<Segment> s;

        Node* ret = new_nodes(1);
        s.push((Segment){0, _size, 1, ret});

        while (!s.empty()) {
            const int begin = s.top().begin;
            const int end = s.top().end;
            const int level = s.top().level;
            Node* node = s.top().node;
            s.pop();

            RT_ASSERT(end - begin >= 2);
            if (end - begin == 2) {
                Node* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1]);
                memcpy(node, _, sizeof(Node));
                delete_nodes(_, 1);
            } else {
                T* keys = _keys + begin;
                P* values = _values + begin;
                const int size = end - begin;
                const int BUILD_GAP_CNT = compute_gap_count(size);

                node->is_two = 0;
                node->build_size = size;
                node->size = size;
                node->fixed = 0;
                node->num_inserts = node->num_insert_to_data = 0;

                int mid1_pos = (size - 1) / 3;
                int mid2_pos = (size - 1) * 2 / 3;

                RT_ASSERT(0 <= mid1_pos);
                RT_ASSERT(mid1_pos < mid2_pos);
                RT_ASSERT(mid2_pos < size - 1);

                const long double mid1_key =
                        (static_cast<long double>(keys[mid1_pos]) + static_cast<long double>(keys[mid1_pos + 1])) / 2;
                const long double mid2_key =
                        (static_cast<long double>(keys[mid2_pos]) + static_cast<long double>(keys[mid2_pos + 1])) / 2;

                node->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

                node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                node->model.b = mid1_target - node->model.a * mid1_key;
                RT_ASSERT(isfinite(node->model.a));
                RT_ASSERT(isfinite(node->model.b));

                const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
                node->model.b += lr_remains;
                node->num_items += lr_remains * 2;

                if (size > 1e6) {
                    node->fixed = 1;
                }

                node->items = new_items(node->num_items);
                const int bitmap_size = BITMAP_SIZE(node->num_items);
                node->none_bitmap = new_bitmap(bitmap_size);
                node->child_bitmap = new_bitmap(bitmap_size);
                memset(node->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
                memset(node->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

                //memset(node->poi_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

                for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {
                    int next = offset + 1, next_i = -1;
                    while (next < size) {
                        next_i = PREDICT_POS(node, keys[next]);
                        if (next_i == item_i) {
                            next ++;
                        } else {
                            break;
                        }
                    }
                    if (next == offset + 1) {
                        BITMAP_CLEAR(node->none_bitmap, item_i);
                        node->items[item_i].comp.data.key = keys[offset];
                        node->items[item_i].comp.data.value = values[offset];
                    } else {
                        // ASSERT(next - offset <= (size+2) / 3);
                        BITMAP_CLEAR(node->none_bitmap, item_i);
                        BITMAP_SET(node->child_bitmap, item_i);
                        node->items[item_i].comp.child = new_nodes(1);
                        s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child});
                    }
                    if (next >= size) {
                        break;
                    } else {
                        item_i = next_i;
                        offset = next;
                    }
                }
            }
        }

        return ret;
    }
    /// bulk build, _keys must be sorted in asc order.
    /// FMCD method.
    Node* build_tree_bulk_fmcd(T* _keys, P* _values, int _size)
    {
        RT_ASSERT(_size > 1);

        typedef struct {
            int begin;
            int end;
            int level; // top level = 1
            Node* node;
        } Segment;
        std::stack<Segment> s;

        Node* ret = new_nodes(1);
        s.push((Segment){0, _size, 1, ret});

        while (!s.empty()) {
            const int begin = s.top().begin;
            const int end = s.top().end;
            const int level = s.top().level;
            Node* node = s.top().node;
            s.pop();

            RT_ASSERT(end - begin >= 2);
            if (end - begin == 2) {
                Node* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1]);
                memcpy(node, _, sizeof(Node));
                delete_nodes(_, 1);
            } else {
                T* keys = _keys + begin;
                P* values = _values + begin;
                const int size = end - begin;
                const int BUILD_GAP_CNT = compute_gap_count(size);

                node->is_two = 0;
                node->build_size = size;
                node->size = size;
                node->fixed = 0;
                node->num_inserts = node->num_insert_to_data = 0;

                // FMCD method
                // Here the implementation is a little different with Algorithm 1 in our paper.
                // In Algorithm 1, U_T should be (keys[size-1-D] - keys[D]) / (L - 2).
                // But according to the derivation described in our paper, M.A should be less than 1 / U_T.
                // So we added a small number (1e-6) to U_T.
                // In fact, it has only a negligible impact of the performance.
                {
                    const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
                    int i = 0;
                    int D = 1;
                    RT_ASSERT(D <= size-1-D);
                    double Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                                (static_cast<double>(L - 2)) + 1e-6;
                    while (i < size - 1 - D) {
                        while (i + D < size && keys[i + D] - keys[i] >= Ut) {
                            i ++;
                        }
                        if (i + D >= size) {
                            break;
                        }
                        D = D + 1;
                        if (D * 3 > size) break;
                        RT_ASSERT(D <= size-1-D);
                        Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                             (static_cast<double>(L - 2)) + 1e-6;
                    }
                    if (D * 3 <= size) {
                        stats.fmcd_success_times ++;

                        node->model.a = 1.0 / Ut;
                        node->model.b = (L - node->model.a * (static_cast<long double>(keys[size - 1 - D]) +
                                                              static_cast<long double>(keys[D]))) / 2;
                        RT_ASSERT(isfinite(node->model.a));
                        RT_ASSERT(isfinite(node->model.b));
                        node->num_items = L;
                    } else {
                        stats.fmcd_broken_times ++;

                        int mid1_pos = (size - 1) / 3;
                        int mid2_pos = (size - 1) * 2 / 3;

                        RT_ASSERT(0 <= mid1_pos);
                        RT_ASSERT(mid1_pos < mid2_pos);
                        RT_ASSERT(mid2_pos < size - 1);

                        const long double mid1_key = (static_cast<long double>(keys[mid1_pos]) +
                                                      static_cast<long double>(keys[mid1_pos + 1])) / 2;
                        const long double mid2_key = (static_cast<long double>(keys[mid2_pos]) +
                                                      static_cast<long double>(keys[mid2_pos + 1])) / 2;

                        node->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                        const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                        const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

                        node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                        node->model.b = mid1_target - node->model.a * mid1_key;
                        RT_ASSERT(isfinite(node->model.a));
                        RT_ASSERT(isfinite(node->model.b));
                    }
                }
                RT_ASSERT(node->model.a >= 0);
                const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
                node->model.b += lr_remains;
                node->num_items += lr_remains * 2;

                if (size > 1e6) {
                    node->fixed = 1;
                }

                node->items = new_items(node->num_items);
                const int bitmap_size = BITMAP_SIZE(node->num_items);
                node->none_bitmap = new_bitmap(bitmap_size);
                node->child_bitmap = new_bitmap(bitmap_size);
                memset(node->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
                memset(node->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

                node->poi_bitmap = new_bitmap(bitmap_size);
                memset(node->poi_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

                //This is where they set the new ones. 
                for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {
                    int next = offset + 1, next_i = -1;
                    while (next < size) {
                        next_i = PREDICT_POS(node, keys[next]);
                        //If the same location then do this. I could check if this is poisoned then can skip. that and check the next original.
                        //If it is JUST poisoned then do not create a child. If it also has original data THEN create the child.
                        if (next_i == item_i) {
                            next ++;
                        } else {
                            break;
                        }
                    }
                    //If poisoned do not add this only if it is original.
                    if (next == offset + 1) {
                        BITMAP_CLEAR(node->none_bitmap, item_i);
                        node->items[item_i].comp.data.key = keys[offset];
                        node->items[item_i].comp.data.value = values[offset];
                    } else {
                        // ASSERT(next - offset <= (size+2) / 3);
                        BITMAP_CLEAR(node->none_bitmap, item_i);
                        BITMAP_SET(node->child_bitmap, item_i);
                        node->items[item_i].comp.child = new_nodes(1);
                        s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child});
                    }
                    if (next >= size) {
                        break;
                    } else {
                        item_i = next_i;
                        offset = next;
                    }
                }
            }
        }

        return ret;
    }
    //=============================================================

    /// bulk build, _keys must be sorted in asc order.
    /// FMCD method.
    //Without adding the poisoning values
    Node* poi_build_tree_bulk_fmcd(T* _keys, P* _values, int _size, T* _keys_real, int _size_real)
    {
        RT_ASSERT(_size > 1);

        typedef struct {
            int begin;
            int end;
            int level; // top level = 1
            Node* node;
            int idx;
            int poi_count;
        } Segment;
        std::stack<Segment> s;

        Node* ret = new_nodes(1);
        s.push((Segment){0, _size, 1, ret, 0, (_size -_size_real)});
        

        while (!s.empty()) {
            const int begin = s.top().begin;
            const int end = s.top().end;
            const int level = s.top().level;
            Node* node = s.top().node;

            int idx = s.top().idx;
            int poi_count_r = s.top().poi_count;

            s.pop();

            RT_ASSERT(end - begin >= 2);
            if (end - begin == 2) {
                //std::cout<< "build 2" <<std::endl;
                Node* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1]);
                memcpy(node, _, sizeof(Node));
                delete_nodes(_, 1);
            } else {
                //std::cout<< "build all" <<std::endl;
                T* keys = _keys + begin;
                P* values = _values + begin;
                const int size = end - begin;
                const int BUILD_GAP_CNT = compute_gap_count(size);

                node->is_two = 0;
                node->build_size = size;
                node->size = size - poi_count_r;
                //node->size = _size_real;
                node->fixed = 0;
                node->num_inserts = node->num_insert_to_data = 0;

                // FMCD method
                // Here the implementation is a little different with Algorithm 1 in our paper.
                // In Algorithm 1, U_T should be (keys[size-1-D] - keys[D]) / (L - 2).
                // But according to the derivation described in our paper, M.A should be less than 1 / U_T.
                // So we added a small number (1e-6) to U_T.
                // In fact, it has only a negligible impact of the performance.
                {
                    const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
                    int i = 0;
                    int D = 1;
                    RT_ASSERT(D <= size-1-D);
                    double Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                                (static_cast<double>(L - 2)) + 1e-6;
                    while (i < size - 1 - D) {
                        while (i + D < size && keys[i + D] - keys[i] >= Ut) {
                            i ++;
                        }
                        if (i + D >= size) {
                            break;
                        }
                        D = D + 1;
                        if (D * 3 > size) break;
                        RT_ASSERT(D <= size-1-D);
                        Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                             (static_cast<double>(L - 2)) + 1e-6;
                    }
                    if (D * 3 <= size) {
                        stats.fmcd_success_times ++;

                        node->model.a = 1.0 / Ut;
                        node->model.b = (L - node->model.a * (static_cast<long double>(keys[size - 1 - D]) +
                                                              static_cast<long double>(keys[D]))) / 2;
                        RT_ASSERT(isfinite(node->model.a));
                        RT_ASSERT(isfinite(node->model.b));
                        node->num_items = L;
                    } else {
                        stats.fmcd_broken_times ++;

                        int mid1_pos = (size - 1) / 3;
                        int mid2_pos = (size - 1) * 2 / 3;

                        RT_ASSERT(0 <= mid1_pos);
                        RT_ASSERT(mid1_pos < mid2_pos);
                        RT_ASSERT(mid2_pos < size - 1);

                        const long double mid1_key = (static_cast<long double>(keys[mid1_pos]) +
                                                      static_cast<long double>(keys[mid1_pos + 1])) / 2;
                        const long double mid2_key = (static_cast<long double>(keys[mid2_pos]) +
                                                      static_cast<long double>(keys[mid2_pos + 1])) / 2;

                        node->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                        const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                        const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

                        node->model.a = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                        node->model.b = mid1_target - node->model.a * mid1_key;
                        RT_ASSERT(isfinite(node->model.a));
                        RT_ASSERT(isfinite(node->model.b));
                    }
                }
                RT_ASSERT(node->model.a >= 0);
                const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
                node->model.b += lr_remains;
                node->num_items += lr_remains * 2;

                if (size > 1e6) {
                    node->fixed = 1;
                }

                node->items = new_items(node->num_items);
                const int bitmap_size = BITMAP_SIZE(node->num_items);
                node->none_bitmap = new_bitmap(bitmap_size);
                node->child_bitmap = new_bitmap(bitmap_size);

                memset(node->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
                memset(node->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);
                
                node->poi_bitmap = new_bitmap(bitmap_size);
                memset(node->poi_bitmap, 0, sizeof(bitmap_t) * bitmap_size);

                //This is where they set the new ones. 
                
                bool real_next = false;
                bool real = false;
                bool real_prev = false;
                int next_prev = 0;
                int count = 0;
                int poi_count = 0;
                int id = 0;

                // for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {

                //     real_next = false;
                //     real_prev = false;
                //     count = 0;
                //     poi_count = 0;

                //     int next = offset + 1, next_i = -1;

                //     if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                //             real_prev = true;
                //             count++;
                //             id = offset;
                //         }
                //         else{
                //             poi_count++;
                //         }

                //     //int next_real = offset + 1;
                //     //The ISSUE IS WITH NEXT_I
                //     while (next < size) {
                        
                //         next_i = PREDICT_POS(node, keys[next]);


                //         //If the same location then do this. I could check if this is poisoned then can skip. that and check the next original.
                //         //If it is JUST poisoned then do not create a child. If it also has original data THEN create the child.
                //         real = false;
                //         if(std::find(_keys_real, _keys_real + _size_real, keys[next]) != _keys_real + _size_real){
                //                 real = true;   
                //                 //idx++;
                //             }
                //         if (next_i == item_i) {
                //             if(real){
                //                 real_next = true;
                //                 count++;
                //                 id = next;
                //                // next_real++;
                //             }
                //             else{
                //                 poi_count++;
                //             }
                            
                //             // else if(real_prev){
                //             //     count++;
                //             // }
                //             next ++;
                //         } else {
                //             break;
                //         }
                //     }
                //     //If poisoned do not add this only if it is original.
                //     if (next == offset + 1) {
                //          //No conflicts and if it is in the real dataset add
                //         // if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                //         if(real_prev){
                //             BITMAP_CLEAR(node->none_bitmap, item_i);
                //             node->items[item_i].comp.data.key = keys[offset];
                //             node->items[item_i].comp.data.value = values[offset];
                //         }
                        
                //     } else {
                //         // ASSERT(next - offset <= (size+2) / 3);
                //         //In the conficts atleast one is real key. Then move to the next.
                //         //DOUBLE CHECK THIS
                //         if(count > 1){
                //             BITMAP_CLEAR(node->none_bitmap, item_i);
                //             BITMAP_SET(node->child_bitmap, item_i);
                //             node->items[item_i].comp.child = new_nodes(1);
                //             s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child, idx, poi_count});
                //         }
                //         else if(count == 1){
                //             if(std::find(_keys_real, _keys_real + _size_real, keys[id]) != _keys_real + _size_real){
                //                 BITMAP_CLEAR(node->none_bitmap, item_i);
                //                 node->items[item_i].comp.data.key = keys[id];
                //                 node->items[item_i].comp.data.value = values[id];
                //             }
                            
                //         }
                //         //If there are conflicts but all are poisoned values then
                //         else{
                //             // if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                //             if(real_prev){
                //                 BITMAP_CLEAR(node->none_bitmap, item_i);
                //                 node->items[item_i].comp.data.key = keys[offset];
                //                 node->items[item_i].comp.data.value = values[offset];
                //             }
                //         }
                       
                //     }
                //     if (next >= size) {
                //         break;
                //     } else {
                //         item_i = next_i;
                //         offset = next;
                //         next_prev = next;
                //     }
                // }

                //THIS WILL SAVE ALL DATA BUT AS SEPERATE TYPES
                //Save all including poisoned. So if there are crashes then it will save them here too
                
                for (int item_i = PREDICT_POS(node, keys[0]), offset = 0; offset < size; ) {

                    real_next = false;
                    real_prev = false;
                    count = 0;
                    poi_count = 0;

                    int next = offset + 1, next_i = -1;

                    if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                            real_prev = true;
                            count++;
                            id = offset;
                        }
                        else{
                            poi_count++;
                        }

                    //int next_real = offset + 1;
                    //The ISSUE IS WITH NEXT_I
                    while (next < size) {
                        
                        next_i = PREDICT_POS(node, keys[next]);


                        //If the same location then do this. I could check if this is poisoned then can skip. that and check the next original.
                        //If it is JUST poisoned then do not create a child. If it also has original data THEN create the child.
                        real = false;
                        if(std::find(_keys_real, _keys_real + _size_real, keys[next]) != _keys_real + _size_real){
                                real = true;   
                                //idx++;
                            }
                        if (next_i == item_i) {
                            if(real){
                                real_next = true;
                                count++;
                                id = next;
                               // next_real++;
                            }
                            else{
                                poi_count++;
                            }
                            
                            // else if(real_prev){
                            //     count++;
                            // }
                            next ++;
                        } else {
                            break;
                        }
                    }
                    //If poisoned do not add this only if it is original.
                    if (next == offset + 1) {
                         //No conflicts and if it is in the real dataset add
                        // if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                        if(real_prev){
                            BITMAP_CLEAR(node->none_bitmap, item_i);
                            node->items[item_i].comp.data.key = keys[offset];
                            node->items[item_i].comp.data.value = values[offset];
                        }
                        else{
                            //ADD SOMETHING
                            BITMAP_CLEAR(node->none_bitmap, item_i);
                            BITMAP_SET(node->poi_bitmap, item_i);
                            
                            node->items[item_i].comp.data.key = keys[offset];
                            node->items[item_i].comp.data.value = values[offset];
                            node->size++;
                            node->poi_size++;
                        }
                        
                    } else {
                        // ASSERT(next - offset <= (size+2) / 3);
                        //In the conficts atleast one is real key. Then move to the next.
                        //DOUBLE CHECK THIS
                        if(count > 1){
                            BITMAP_CLEAR(node->none_bitmap, item_i);
                            BITMAP_SET(node->child_bitmap, item_i);
                            node->items[item_i].comp.child = new_nodes(1);
                            s.push((Segment){begin + offset, begin + next, level + 1, node->items[item_i].comp.child, idx, poi_count});
                        }
                        else if(count == 1){
                            if(std::find(_keys_real, _keys_real + _size_real, keys[id]) != _keys_real + _size_real){
                                BITMAP_CLEAR(node->none_bitmap, item_i);
                                node->items[item_i].comp.data.key = keys[id];
                                node->items[item_i].comp.data.value = values[id];
                            }
                            
                        }
                        //If there are conflicts but all are poisoned values then
                        else{
                            // if(std::find(_keys_real, _keys_real + _size_real, keys[offset]) != _keys_real + _size_real){
                            if(real_prev){
                                BITMAP_CLEAR(node->none_bitmap, item_i);
                                node->items[item_i].comp.data.key = keys[offset];
                                node->items[item_i].comp.data.value = values[offset];
                            }
                        }
                       
                    }
                    if (next >= size) {
                        break;
                    } else {
                        item_i = next_i;
                        offset = next;
                        next_prev = next;
                    }
                }
            }
        }

        return ret;
    }

    void destory_pending()
    {
        while (!pending_two.empty()) {
            Node* node = pending_two.top(); pending_two.pop();

            delete_items(node->items, node->num_items);
            const int bitmap_size = BITMAP_SIZE(node->num_items);
            delete_bitmap(node->none_bitmap, bitmap_size);
            delete_bitmap(node->child_bitmap, bitmap_size);

            delete_bitmap(node->poi_bitmap, bitmap_size);
            delete_nodes(node, 1);
        }
    }

    void destroy_tree(Node* root)
    {
        std::stack<Node*> s;
        s.push(root);
        while (!s.empty()) {
            Node* node = s.top(); s.pop();

            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->child_bitmap, i) == 1) {
                    s.push(node->items[i].comp.child);
                }
            }

            if (node->is_two) {
                RT_ASSERT(node->build_size == 2);
                RT_ASSERT(node->num_items == 8);
                node->size = 2;
                node->num_inserts = node->num_insert_to_data = 0;
                node->none_bitmap[0] = 0xff;
                node->child_bitmap[0] = 0;
                pending_two.push(node);

                node->poi_bitmap[0] = 0;
            } else {
                delete_items(node->items, node->num_items);
                const int bitmap_size = BITMAP_SIZE(node->num_items);
                delete_bitmap(node->none_bitmap, bitmap_size);
                delete_bitmap(node->child_bitmap, bitmap_size);

                delete_bitmap(node->poi_bitmap, bitmap_size);
                delete_nodes(node, 1);
            }
        }
    }
    //Useful to get the keys and values for the subtree of a node
    void scan_and_destory_tree(Node* _root, T* keys, P* values, bool destory = true)
    {
        //Save the real keys in a seperate vector. Use node-size and node-poi_size to determine the size
        typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, _root));
        while (!s.empty()) {
            int begin = s.top().first;
            Node* node = s.top().second;
            const int SHOULD_END_POS = begin + node->size;
            s.pop();

            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (BITMAP_GET(node->child_bitmap, i) == 0) {
                        keys[begin] = node->items[i].comp.data.key;
                        values[begin] = node->items[i].comp.data.value;
                        begin ++;
                    } else {
                        //std::cout << " 1 " << begin << std::endl;
                        s.push(Segment(begin, node->items[i].comp.child));
                        begin += node->items[i].comp.child->size;
                        //std::cout << " 2 " << begin << std::endl;
                    }
                }
            }
            //std::cout << SHOULD_END_POS << " , " << begin << std::endl;
            RT_ASSERT(SHOULD_END_POS == begin);
            

            if (destory) {
                if (node->is_two) {
                    RT_ASSERT(node->build_size == 2);
                    RT_ASSERT(node->num_items == 8);
                    node->size = 2;
                    node->num_inserts = node->num_insert_to_data = 0;
                    node->none_bitmap[0] = 0xff;
                    node->child_bitmap[0] = 0;

                    node->poi_bitmap[0] = 0;
                    pending_two.push(node);
                } else {
                    delete_items(node->items, node->num_items);
                    const int bitmap_size = BITMAP_SIZE(node->num_items);
                    delete_bitmap(node->none_bitmap, bitmap_size);
                    delete_bitmap(node->child_bitmap, bitmap_size);

                    delete_bitmap(node->poi_bitmap, bitmap_size);
                    delete_nodes(node, 1);
                }
            }
        }
    }

    void scan_and_destory_tree2(Node* _root, T* keys, P* values, T* keys_real, bool destory = true)
    {
        //Save the real keys in a seperate vector. Use node-size and node-poi_size to determine the size
        typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, _root));
        int p = 0;

        while (!s.empty()) {
            int begin = s.top().first;
            Node* node = s.top().second;
            const int SHOULD_END_POS = begin + node->size;
            s.pop();

            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (BITMAP_GET(node->child_bitmap, i) == 0) {
                        keys[begin] = node->items[i].comp.data.key;
                        values[begin] = node->items[i].comp.data.value;
                        begin ++;

                        if(BITMAP_GET(node->poi_bitmap, i) == 0){
                            keys_real[p] = node->items[i].comp.data.key;
                            p++;
                        }
                    } else {
                        //std::cout << " 1 " << begin << std::endl;
                        s.push(Segment(begin, node->items[i].comp.child));
                        begin += node->items[i].comp.child->size;
                        //std::cout << " 2 " << begin << std::endl;
                    }
                }
            }
            //std::cout << SHOULD_END_POS << " , " << begin << std::endl;
            RT_ASSERT(SHOULD_END_POS == begin);
            

            if (destory) {
                if (node->is_two) {
                    RT_ASSERT(node->build_size == 2);
                    RT_ASSERT(node->num_items == 8);
                    node->size = 2;
                    node->num_inserts = node->num_insert_to_data = 0;
                    node->none_bitmap[0] = 0xff;
                    node->child_bitmap[0] = 0;

                    node->poi_bitmap[0] = 0;
                    pending_two.push(node);
                } else {
                    delete_items(node->items, node->num_items);
                    const int bitmap_size = BITMAP_SIZE(node->num_items);
                    delete_bitmap(node->none_bitmap, bitmap_size);
                    delete_bitmap(node->child_bitmap, bitmap_size);

                    delete_bitmap(node->poi_bitmap, bitmap_size);
                    delete_nodes(node, 1);
                }
            }
        }
    }

    void scan_subtree(Node* _root, std::vector<T> * keys, std::vector<P> * values, bool destory = false)
    {
        typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, _root));
        while (!s.empty()) {
            int begin = s.top().first;
            Node* node = s.top().second;
            const int SHOULD_END_POS = begin + node->size;
            s.pop();

            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (BITMAP_GET(node->child_bitmap, i) == 0) {
                        keys->push_back(node->items[i].comp.data.key);
                        values->push_back(node->items[i].comp.data.value);
                        begin ++;
                    } else {
                        s.push(Segment(begin, node->items[i].comp.child));
                        begin += node->items[i].comp.child->size;
                    }
                }
            }
            RT_ASSERT(SHOULD_END_POS == begin);

            if (destory) {
                if (node->is_two) {
                    RT_ASSERT(node->build_size == 2);
                    RT_ASSERT(node->num_items == 8);
                    node->size = 2;
                    node->num_inserts = node->num_insert_to_data = 0;
                    node->none_bitmap[0] = 0xff;
                    node->child_bitmap[0] = 0;

                    node->poi_bitmap[0] = 0;
                    pending_two.push(node);
                } else {
                    delete_items(node->items, node->num_items);
                    const int bitmap_size = BITMAP_SIZE(node->num_items);
                    delete_bitmap(node->none_bitmap, bitmap_size);
                    delete_bitmap(node->child_bitmap, bitmap_size);

                    delete_bitmap(node->poi_bitmap, bitmap_size);
                    delete_nodes(node, 1);
                }
            }
        }
    }

    void scan_subtree2(Node* _root, T* keys, P* values, bool destory = false)
    {
        typedef std::pair<int, Node*> Segment; // <begin, Node*>
        std::stack<Segment> s;

        s.push(Segment(0, _root));
        while (!s.empty()) {
            int begin = s.top().first;
            Node* node = s.top().second;
            const int SHOULD_END_POS = begin + node->size;
            s.pop();
            //std::cout << "Here scan1 " << begin <<" " << SHOULD_END_POS<< std::endl;
            for (int i = 0; i < node->num_items; i ++) {
                if (BITMAP_GET(node->none_bitmap, i) == 0) {
                    if (BITMAP_GET(node->child_bitmap, i) == 0) {
                        keys[begin] = node->items[i].comp.data.key;
                        values[begin] = node->items[i].comp.data.value;
                        begin ++;
                    } else if(BITMAP_GET(node->child_bitmap, i) == 1){
                        s.push(Segment(begin, node->items[i].comp.child));
                        begin += node->items[i].comp.child->size;
                        //std::cout << "Here scan2 " << begin <<" " << SHOULD_END_POS<< std::endl;
                    }
                }
            }
            //std::cout << "Here scan " << begin <<" " << SHOULD_END_POS<< std::endl;
            RT_ASSERT(SHOULD_END_POS == begin);

            if (destory) {
                if (node->is_two) {
                    RT_ASSERT(node->build_size == 2);
                    RT_ASSERT(node->num_items == 8);
                    node->size = 2;
                    node->num_inserts = node->num_insert_to_data = 0;
                    node->none_bitmap[0] = 0xff;
                    node->child_bitmap[0] = 0;
                    pending_two.push(node);
                } else {
                    delete_items(node->items, node->num_items);
                    const int bitmap_size = BITMAP_SIZE(node->num_items);
                    delete_bitmap(node->none_bitmap, bitmap_size);
                    delete_bitmap(node->child_bitmap, bitmap_size);
                    delete_nodes(node, 1);
                }
            }
        }
    }

    //Has the adjustment here
    //here need to check if the poisoned or not before entering. Basically use the poi_bulk_load_tree here itself
    Node* insert_tree(Node* _node, const T& key, const P& value)
    {
        constexpr int MAX_DEPTH = 128;
        Node* path[MAX_DEPTH];
        int path_size = 0;
        int insert_to_data = 0;

        for (Node* node = _node; ; ) {
            RT_ASSERT(path_size < MAX_DEPTH);
            path[path_size ++] = node;

            node->size ++;
            node->num_inserts ++;
            int pos = PREDICT_POS(node, key);
            if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                BITMAP_CLEAR(node->none_bitmap, pos);
                node->items[pos].comp.data.key = key;
                node->items[pos].comp.data.value = value;
                break;
            } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
                BITMAP_SET(node->child_bitmap, pos);
                node->items[pos].comp.child = build_tree_two(key, value, node->items[pos].comp.data.key, node->items[pos].comp.data.value);
                insert_to_data = 1;
                break;
            } else {
                node = node->items[pos].comp.child;
            }
        }
        for (int i = 0; i < path_size; i ++) {
            path[i]->num_insert_to_data += insert_to_data;
        }

        for (int i = 0; i < path_size; i ++) {
            Node* node = path[i];
            const int num_inserts = node->num_inserts;
            const int num_insert_to_data = node->num_insert_to_data;
            const bool need_rebuild = node->fixed == 0 && node->size >= node->build_size * 4 && node->size >= 64 && num_insert_to_data * 10 >= num_inserts;
            
            if (need_rebuild) {
                std::cout << "rebuild" << std::endl;
                const int ESIZE = node->size;
                T* keys = new T[ESIZE];
                P* values = new P[ESIZE];

                #if COLLECT_TIME
                auto start_time_scan = std::chrono::high_resolution_clock::now();
                #endif
                scan_and_destory_tree(node, keys, values);
                #if COLLECT_TIME
                auto end_time_scan = std::chrono::high_resolution_clock::now();
                auto duration_scan = end_time_scan - start_time_scan;
                stats.time_scan_and_destory_tree += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_scan).count() * 1e-9;
                #endif

                #if COLLECT_TIME
                auto start_time_build = std::chrono::high_resolution_clock::now();
                #endif
                Node* new_node = build_tree_bulk(keys, values, ESIZE);
                #if COLLECT_TIME
                auto end_time_build = std::chrono::high_resolution_clock::now();
                auto duration_build = end_time_build - start_time_build;
                stats.time_build_tree_bulk += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_build).count() * 1e-9;
                #endif

                delete[] keys;
                delete[] values;

                path[i] = new_node;
                if (i > 0) {
                    int pos = PREDICT_POS(path[i-1], key);
                    path[i-1]->items[pos].comp.child = new_node;
                }

                break;
            }
        }

        return path[0];
    }

    Node* insert_tree2(Node* _node, const T& key, const P& value)
    {
        constexpr int MAX_DEPTH = 128;
        Node* path[MAX_DEPTH];
        int path_size = 0;
        int insert_to_data = 0;

        for (Node* node = _node; ; ) {
            RT_ASSERT(path_size < MAX_DEPTH);
            path[path_size ++] = node;

            node->size ++;
            node->num_inserts ++;
            int pos = PREDICT_POS(node, key);
            if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                BITMAP_CLEAR(node->none_bitmap, pos);
                node->items[pos].comp.data.key = key;
                node->items[pos].comp.data.value = value;
                break;
            } 
            else if(BITMAP_GET(node->poi_bitmap, pos) == 1){
                BITMAP_CLEAR(node->poi_bitmap, pos);
                node->items[pos].comp.data.key = key;
                node->items[pos].comp.data.value = value;
                break;
            }
            
            else if (BITMAP_GET(node->child_bitmap, pos) == 0 && BITMAP_GET(node->poi_bitmap, pos) == 0) {
                BITMAP_SET(node->child_bitmap, pos);
                node->items[pos].comp.child = build_tree_two(key, value, node->items[pos].comp.data.key, node->items[pos].comp.data.value);
                insert_to_data = 1;
                break;
            } else {
                node = node->items[pos].comp.child;
            }
        }
        for (int i = 0; i < path_size; i ++) {
            path[i]->num_insert_to_data += insert_to_data;
        }

        for (int i = 0; i < path_size; i ++) {
            Node* node = path[i];
            const int num_inserts = node->num_inserts;
            const int num_insert_to_data = node->num_insert_to_data;
            const bool need_rebuild = node->fixed == 0 && node->size >= node->build_size * 4 && node->size >= 64 && num_insert_to_data * 10 >= num_inserts ;
            
            if (need_rebuild) {
                //std::cout << "rebuild" << std::endl;
                const int ESIZE = node->size;
                T* keys = new T[ESIZE];
                P* values = new P[ESIZE];

                T* keys_real = new T[ESIZE - node->poi_size];

                #if COLLECT_TIME
                auto start_time_scan = std::chrono::high_resolution_clock::now();
                #endif
                scan_and_destory_tree2(node, keys, values, keys_real);
                #if COLLECT_TIME
                auto end_time_scan = std::chrono::high_resolution_clock::now();
                auto duration_scan = end_time_scan - start_time_scan;
                stats.time_scan_and_destory_tree += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_scan).count() * 1e-9;
                #endif

                #if COLLECT_TIME
                auto start_time_build = std::chrono::high_resolution_clock::now();
                #endif
                // Node* new_node = build_tree_bulk(keys, values, ESIZE);
                Node* new_node = poi_build_tree_bulk(keys, values, ESIZE,keys_real, ESIZE-node->poi_size);
                #if COLLECT_TIME
                auto end_time_build = std::chrono::high_resolution_clock::now();
                auto duration_build = end_time_build - start_time_build;
                stats.time_build_tree_bulk += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_build).count() * 1e-9;
                #endif

                delete[] keys;
                delete[] values;

                path[i] = new_node;
                if (i > 0) {
                    int pos = PREDICT_POS(path[i-1], key);
                    path[i-1]->items[pos].comp.child = new_node;
                }

                break;
            }
        }

        return path[0];
    }

    //To create new ones with the inverse poisoned values
    Node* poison_bulk_load(Node* _node, Node* node_poi, const T& key, const P& value, T* keys,P* values, int ESIZE_POI)
    {
        //  constexpr int MAX_DEPTH = 128;
        // Node* path[MAX_DEPTH];
        // int path_size = 0;
        // int insert_to_data = 0;

        // for (Node* node = _node; ; ) {
        //     RT_ASSERT(path_size < MAX_DEPTH);
        //     path[path_size ++] = node;

        //     node->size ++;
        //     node->num_inserts ++;
        //     int pos = PREDICT_POS(node, key);
        //     if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        //         //BITMAP_CLEAR(node->none_bitmap, pos);
        //         //node->items[pos].comp.data.key = key;
        //         //node->items[pos].comp.data.value = value;
        //         break;
        //     } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        //         //BITMAP_SET(node->child_bitmap, pos);
        //         //node->items[pos].comp.child = build_tree_two(key, value, node->items[pos].comp.data.key, node->items[pos].comp.data.value);
        //         //insert_to_data = 1;
        //         break;
        //     } else {
        //         node = node->items[pos].comp.child;
        //     }
        // }
        constexpr int MAX_DEPTH = 128;
        Node* path[MAX_DEPTH];
        int path_size = 0;
        Node* node = _node;
        // path[path_size] = node;
        // path_size++;

        while (true) {
            
            int pos = PREDICT_POS(node, key);
            path[path_size] = node;
            //std::cout<< path[path_size] <<" , ";
            //std::cout<<"=========== " << path_size << std::endl;

            path_size++;

            if (BITMAP_GET(node->child_bitmap, pos) == 1) {
                node = node->items[pos].comp.child;
                
            } else {
                if (false) {
                    break;
                    //return node->items[pos].comp.data.value;
                } else {
                    if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                        RT_ASSERT(false);
                    } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
                        RT_ASSERT(node->items[pos].comp.data.key == key);
                        break;
                        //return node->items[pos].comp.data.value;
                    }
                }
            }
        }

        // for (int i = 0; i < path_size; i ++) {
        //     //path[i]->num_insert_to_data += insert_to_data;
        //     std::cout<< path[i] <<" , ";
        // }

        for (int i = 0; i < 1; i ++) {
            Node* node = path[path_size-1];

            if(node != node_poi){
                std::cout<< "NOT SAME" <<std::endl;
            }
            // std::cout<< node <<std::endl;
            // std::cout<< node_poi <<std::endl;
            //Node* node = node_poi;
            // const int num_inserts = node->num_inserts;
            // const int num_insert_to_data = node->num_insert_to_data;
            // const bool need_rebuild = node->fixed == 0 && node->size >= node->build_size * 4 && node->size >= 64 && num_insert_to_data * 10 >= num_inserts;

            if (true) {
                const int ESIZE = node->size;
                T* keys2 = new T[ESIZE];
                P* values2 = new P[ESIZE];

                #if COLLECT_TIME
                auto start_time_scan = std::chrono::high_resolution_clock::now();
                #endif
                scan_and_destory_tree(node, keys2, values2);
                #if COLLECT_TIME
                auto end_time_scan = std::chrono::high_resolution_clock::now();
                auto duration_scan = end_time_scan - start_time_scan;
                stats.time_scan_and_destory_tree += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_scan).count() * 1e-9;
                #endif

                #if COLLECT_TIME
                auto start_time_build = std::chrono::high_resolution_clock::now();
                #endif
                Node* new_node = build_tree_bulk(keys, values, ESIZE_POI);
                #if COLLECT_TIME
                auto end_time_build = std::chrono::high_resolution_clock::now();
                auto duration_build = end_time_build - start_time_build;
                stats.time_build_tree_bulk += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_build).count() * 1e-9;
                #endif

                delete[] keys;
                delete[] values;
                //new_node->size = ESIZE_POI;

                //path[i] = new_node;

                //The answer is somewhere here the assigning stage
                // if (i > 0) {
                //     int pos = PREDICT_POS(path[i-1], key);
                //     path[i-1]->items[pos].comp.child = new_node;
                // }
                int pos = PREDICT_POS(path[path_size-2], key);
                path[path_size-2]->items[pos].comp.child = new_node;

                //Increasing the height of the above ones by this amount.
                // path[path_size-2]->size = path[path_size-2]->size+ESIZE_POI-ESIZE;
                //Increasing the size of all the ones above the path.
                for(int j = path_size-2; j > 0 ; j--){
                    path[j]->size = path[j]->size + ESIZE_POI-ESIZE;
                    //path[j]->size = path[j]->size + ESIZE_POI;
                }
                

                break;
            }
        }
        //std::cout<< "Path 0 "<< path[0] <<std::endl;

        return path[0];
    }

    Node* poison_bulk_load_real(Node* _node, Node* node_poi, const T& key, const P& value, T* keys,P* values, int ESIZE_POI, T* keys_real, int _size_real)
    {
        //  constexpr int MAX_DEPTH = 128;
        // Node* path[MAX_DEPTH];
        // int path_size = 0;
        // int insert_to_data = 0;

        // for (Node* node = _node; ; ) {
        //     RT_ASSERT(path_size < MAX_DEPTH);
        //     path[path_size ++] = node;

        //     node->size ++;
        //     node->num_inserts ++;
        //     int pos = PREDICT_POS(node, key);
        //     if (BITMAP_GET(node->none_bitmap, pos) == 1) {
        //         //BITMAP_CLEAR(node->none_bitmap, pos);
        //         //node->items[pos].comp.data.key = key;
        //         //node->items[pos].comp.data.value = value;
        //         break;
        //     } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
        //         //BITMAP_SET(node->child_bitmap, pos);
        //         //node->items[pos].comp.child = build_tree_two(key, value, node->items[pos].comp.data.key, node->items[pos].comp.data.value);
        //         //insert_to_data = 1;
        //         break;
        //     } else {
        //         node = node->items[pos].comp.child;
        //     }
        // }
        constexpr int MAX_DEPTH = 128;
        Node* path[MAX_DEPTH];
        int path_size = 0;
        Node* node = _node;

        //COMMENT IF ISSUE
        Node* new_node;
        // path[path_size] = node;
        // path_size++;

        while (true) {
            
            int pos = PREDICT_POS(node, key);
            path[path_size] = node;
            //std::cout<< path[path_size] <<" , ";
            //std::cout<<"=========== " << path_size << std::endl;

            path_size++;

            if (BITMAP_GET(node->child_bitmap, pos) == 1) {
                node = node->items[pos].comp.child;
                
            } else {
                if (false) {
                    break;
                    //return node->items[pos].comp.data.value;
                } else {
                    if (BITMAP_GET(node->none_bitmap, pos) == 1) {
                        RT_ASSERT(false);
                    } else if (BITMAP_GET(node->child_bitmap, pos) == 0) {
                        RT_ASSERT(node->items[pos].comp.data.key == key);
                        break;
                        //return node->items[pos].comp.data.value;
                    }
                }
            }
        }

        // for (int i = 0; i < path_size; i ++) {
        //     //path[i]->num_insert_to_data += insert_to_data;
        //     std::cout<< path[i] <<" , ";
        // }

        for (int i = 0; i < 1; i ++) {
            Node* node = path[path_size-1];

            if(node != node_poi){
                std::cout<< "NOT SAME" <<std::endl;
            }
            // std::cout<< node <<std::endl;
            // std::cout<< node_poi <<std::endl;
            //Node* node = node_poi;
            // const int num_inserts = node->num_inserts;
            // const int num_insert_to_data = node->num_insert_to_data;
            // const bool need_rebuild = node->fixed == 0 && node->size >= node->build_size * 4 && node->size >= 64 && num_insert_to_data * 10 >= num_inserts;
            
            if (true) {
                const int ESIZE = node->size;
                T* keys2 = new T[ESIZE];
                P* values2 = new P[ESIZE];

                #if COLLECT_TIME
                auto start_time_scan = std::chrono::high_resolution_clock::now();
                #endif
                scan_and_destory_tree(node, keys2, values2);
                #if COLLECT_TIME
                auto end_time_scan = std::chrono::high_resolution_clock::now();
                auto duration_scan = end_time_scan - start_time_scan;
                stats.time_scan_and_destory_tree += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_scan).count() * 1e-9;
                #endif

                #if COLLECT_TIME
                auto start_time_build = std::chrono::high_resolution_clock::now();
                #endif
                //UNCOMMENT IF ISSUE
                //Node* new_node = poi_build_tree_bulk(keys, values, ESIZE_POI,keys_real, _size_real);
                //std::cout<< " SAME 1 " <<std::endl;
                new_node = poi_build_tree_bulk(keys, values, ESIZE_POI,keys_real, _size_real);
                //std::cout<< " SAME 2 " <<std::endl;
                #if COLLECT_TIME
                auto end_time_build = std::chrono::high_resolution_clock::now();
                auto duration_build = end_time_build - start_time_build;
                stats.time_build_tree_bulk += std::chrono::duration_cast<std::chrono::nanoseconds>(duration_build).count() * 1e-9;
                #endif

                delete[] keys;
                delete[] values;
                //new_node->size = ESIZE_POI;

                //path[i] = new_node;

                //The answer is somewhere here the assigning stage
                // if (i > 0) {
                //     int pos = PREDICT_POS(path[i-1], key);
                //     path[i-1]->items[pos].comp.child = new_node;
                // }
                int pos = PREDICT_POS(path[path_size-2], key);
                path[path_size-2]->items[pos].comp.child = new_node;

                //Increasing the height of the above ones by this amount.
                // path[path_size-2]->size = path[path_size-2]->size+ESIZE_POI-ESIZE;
                //Increasing the size of all the ones above the path.
                for(int j = path_size-2; j > 0 ; j--){
                    //path[j]->size = path[j]->size + ESIZE_POI - ESIZE;
                }
                break;
            }
        }
        //std::cout<< "Path 0 "<< path[0] <<std::endl;
        //UNCOMMENT IF ISSUE
        // return path[0];
        
        return new_node;
    }
};

#endif // __LIPP_H__
