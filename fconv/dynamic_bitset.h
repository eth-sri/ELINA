/*
 * Disclaimer: most of the design and implementation of dynamic_bitset is guided
 * by cddlib's and Boost's implementation of bitset.
 */

#include <limits.h>
#include <vector>
#include "asrt.h"
#include "setoper.h"

using namespace std;

// It is important that block_t is unsigned because in that case
// bit shifting operations work correctly.
using block_t = unsigned long;
using set_t = block_t*;

set_t set_create(int num_bits);

set_t set_resize(set_t set, int num_bits);

void set_free(set_t set);

set_t set_copy(set_t set);

bool set_test_bit(set_t set, int n);

void set_enable_bit(set_t set, int n);

void set_enable_all(set_t set);

bool set_intersect_by_any(set_t first, set_t second);

set_t set_intersect(set_t first, set_t second);

bool set_equal(set_t first, set_t second);

bool set_is_subset_of(set_t potential_subset, set_t potential_superset);

int set_count(set_t set);

int set_size(set_t set);

set_t set_from_cdd(set_type cdd_set);

vector<set_t> set_arr_create(int num_sets, int num_bits);

void set_arr_free(const vector<set_t>& sets);

vector<set_t> set_arr_transpose(const vector<set_t>& input);
