/*
 * Disclaimer: most of the design and implementation of dynamic_bitset is guided
 * by cddlib's and Boost's implementation of bitset.
 */

#include <limits.h>
#include <vector>
#include "asrt.h"

using namespace std;

using block_t = unsigned long;
using set_t = block_t*;

set_t set_create(int num_bits);

void set_free(set_t set);

bool set_test(set_t set, int n);

void set_set(set_t set, int n);

bool set_intersect_by_any(set_t first, set_t second);

set_t set_intersect(set_t first, set_t second);

bool set_is_subset_of(set_t potential_subset, set_t potential_superset);

int set_count(set_t set);

int set_size(set_t set);

vector<set_t> set_transpose(const vector<set_t>& input);

set_t set_resize(set_t set, int n);

// TODO: Implement a function for concatenation of two sets.
