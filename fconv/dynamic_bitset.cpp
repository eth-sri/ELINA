/*
 * Disclaimer: most of the design and implementation of dynamic_bitset is guided
 * by cddlib's and Boost's implementation of bitset.
 */

#include <string.h>
#include <iostream>
#include <cassert>
#include "dynamic_bitset.h"
#include "setoper.h"

using namespace std;

constexpr block_t ONE = 1;
constexpr block_t ALL_SET = -1;
constexpr int BITS_IN_BLOCK = sizeof(block_t) * CHAR_BIT;

constexpr char BYTE2COUNT[256] = {
        0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
        3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

inline int set_number_of_blocks(const int num_bits) {
    assert(num_bits >= 0 && "Expected a non-negative num_bits.");
    return num_bits == 0 ? 1 : (num_bits - 1) / BITS_IN_BLOCK + 2;
}

set_t set_create(const int num_bits)
{
    assert(num_bits >= 0 && "Expected a non-negative num_bits.");
    int num_blocks = set_number_of_blocks(num_bits);
    set_t set = (block_t*) calloc(num_blocks, sizeof(block_t));
    set[0] = (block_t) num_bits;
    for (int i = 1; i < num_blocks; i++) {
        set[i] = 0;
    }
    return set;
}

vector<set_t> set_create_vector(const int num_sets, const int num_bits) {
    assert(num_sets >= 0 && num_bits >= 0 && "Expected a non-negative num_sets and num_bits.");
    vector<set_t> sets(num_sets);
    for (int i = 0; i < num_sets; i++) {
        sets[i] = set_create(num_bits);
    }
    return sets;
}

void set_free(const set_t set)
{
    free(set);
}

void set_free_vector(const vector<set_t>& sets) {
    for (auto set : sets) {
        free(set);
    }
}

set_t set_copy(const set_t set) {
    int num_blocks = set_number_of_blocks(set[0]);
    set_t set_copy = (block_t*) calloc(num_blocks, sizeof(block_t));
    memcpy(set_copy, set, num_blocks * sizeof(block_t));
    return set_copy;
}

bool set_test(const set_t set, const int n)
{
    assert(0 <= n && n < (int) set[0] && "Tested element should be within range of set.");
    int block_i = n / BITS_IN_BLOCK + 1;
    int bit_i = n % BITS_IN_BLOCK;
    block_t block = set[block_i];
    return block & (ONE << bit_i);
}

void set_set(const set_t set, const int n) {
    assert(0 <= n && n < (int) set[0] && "Tested element should be within range of set.");
    int block_i = n / BITS_IN_BLOCK + 1;
    int bit_i = n % BITS_IN_BLOCK;
    set[block_i] = set[block_i] | (ONE << bit_i);
}

void set_set_all(const set_t set) {
    int num_bits = (int) set[0];
    int num_blocks = set_number_of_blocks(num_bits);
    // Note that the last block cannot be set - otherwise the count would not be correct.
    for (int i = 1; i < num_blocks - 1; i++) {
        set[i] = ALL_SET;
    }
    int bits_left = num_bits - BITS_IN_BLOCK * (num_blocks - 2);
    assert(bits_left > 0 && "At least one should be left to set, otherwise number of blocks is incorrect.");

    // This operation sets the bits_left only.
    set[num_blocks - 1] = ALL_SET >> (BITS_IN_BLOCK - bits_left);
}

bool set_intersect_by_any(const set_t first, const set_t second)
{
    assert(first[0] == second[0] && "Sets expected to be of the same size.");

    int num_blocks = set_number_of_blocks(first[0]);
    for (int i = 1; i < num_blocks; i++) {
        if (first[i] & second[i]) {
            return true;
        }
    }
    return false;
}

set_t set_intersect(const set_t first, const set_t second) {
    assert(first[0] == second[0] && "Sets expected to be of the same size.");

    int num_blocks = set_number_of_blocks(first[0]);
    set_t res = (block_t*) calloc(num_blocks, sizeof(block_t));
    res[0] = first[0];

    for (int i = 1; i < num_blocks; i++) {
        res[i] = first[i] & second[i];
    }

    return res;
}

bool set_equal(const set_t first, const set_t second) {
    assert(first[0] == second[0] && "Sets expected to be of the same size.");

    int num_blocks = set_number_of_blocks(first[0]);
    for (int i = 1; i < num_blocks; i++) {
        if (first[i] != second[i]) {
            return false;
        }
    }
    return true;
}

bool set_is_subset_of(const set_t potential_subset, const set_t potential_superset)
{
    assert(potential_superset[0] == potential_subset[0] && "Sets expected to be of the same size.");

    int num_blocks = set_number_of_blocks(potential_superset[0]);
    for (int i = 1; i < num_blocks; i++) {
        if ((potential_superset[i] | potential_subset[i]) != potential_superset[i]) {
            return false;
        }
    }

    return true;
}

int set_count(const set_t set)
{
    unsigned char* byte_p = (unsigned char*) (&set[1]);
    int num_bytes = (set_number_of_blocks(set[0]) - 1) * sizeof(block_t);
    int count = 0;
    for (int i = 0; i < num_bytes; i++) {
        count += BYTE2COUNT[*byte_p];
        byte_p++;
    }
    return count;
}

int set_size(const set_t set)
{
    return (int) set[0];
}

// Maybe has potential for optimization.
vector<set_t> set_transpose(const vector<set_t>& input) {
    if (input.empty()) {
        return {};
    }

    const int rows = (int) input.size();
    const int cols = (int) input[0][0];

    vector<set_t> output(cols);
    for (int i = 0; i < cols; i++) {
        output[i] = set_create(rows);
    }

    for (int i = 0; i < rows; i++) {
        const set_t in_row = input[i];
        assert((int) in_row[0] == cols && "All sets should be of the same size.");
        for (int j = 0; j < cols; j++) {
            if (!set_test(in_row, j)) {
                continue;
            }
            set_set(output[j], i);
        }
    }

    return output;
}

set_t set_resize(set_t set, const int num_bits) {
    const int old_num_bits = set[0];
    if (num_bits == old_num_bits) {
        return set;
    }
    // Update number of bits.
    set[0] = num_bits;
    const int old_blocks = set_number_of_blocks(old_num_bits);
    const int blocks = set_number_of_blocks(num_bits);
    if (old_blocks == blocks) {
        // Same amount of memory allocated - no need to reallocate.
        return set;
    }
    set = (set_t) realloc(set, blocks * sizeof(block_t));
    for (int i = old_blocks; i < blocks; i++) {
        // If any new blocks were allocated - set them to zero.
        set[i] = 0;
    }
    return set;
}

set_t set_from_cdd(set_type cdd_set) {
    // For some reason CDD has one more bit than necessary.
    int cdd_num_bits = cdd_set[0];
    cdd_set[0]--;

    int cdd_blocks = set_blocks(cdd_num_bits);
    // Blocks can only be smaller than cdd_blocks.
    int blocks = set_number_of_blocks(cdd_num_bits - 1);

    if (blocks == cdd_blocks) {
        return (set_t) cdd_set;
    }

    assert(blocks == cdd_blocks - 1 && "blocks should equals cdd_blocks - 1");
    return (set_t) realloc(cdd_set, blocks * sizeof(block_t));
}
