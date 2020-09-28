/*
 * Disclaimer: most of the design and implementation of dynamic_bitset is guided
 * by cddlib's and Boost's implementation of bitset.
 */

#include "dynamic_bitset.h"

using namespace std;

constexpr block_t ONE = 1;
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

void set_free(const set_t set)
{
    free(set);
}

bool set_test(const set_t set, const int n)
{
    assert(0 <= n && n < set[0] && "Tested element should be within range of set.");
    int block_i = n / BITS_IN_BLOCK + 1;
    int bit_i = n % BITS_IN_BLOCK;
    block_t block = set[block_i];
    return block & (ONE << bit_i);
}

void set_set(const set_t set, const int n) {
    assert(0 <= n && n < set[0] && "Tested element should be within range of set.");
    int block_i = n / BITS_IN_BLOCK + 1;
    int bit_i = n % BITS_IN_BLOCK;
    set[block_i] = set[block_i] | (ONE << bit_i);
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
        assert(in_row[0] == cols && "All sets should be of the same size.");
        for (int j = 0; j < cols; j++) {
            if (!set_test(in_row, j)) {
                continue;
            }
            set_set(output[j], i);
        }
    }

    return output;
}

set_t set_resize(set_t set, const int n) {
    const int old_n = set[0];
    if (n == old_n) {
        return set;
    }
    // Update number of bits.
    set[0] = n;
    const int old_blocks = set_number_of_blocks(old_n);
    const int blocks = set_number_of_blocks(old_n);
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
