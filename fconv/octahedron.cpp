#include <set>
#include <map>
#include <cassert>
#include <algorithm>
#include "octahedron.h"
#include "utils.h"
#include "mpq.h"
#include "dynamic_bitset.h"
#include "fp_mat.h"

// The absolutely proper way would be to set it based on the size of the input constraints.
// But I will just keep it fixed for now - I use this value to do partially solve problem in fp precision.
constexpr double DOUBLE_PRECISION = 0.1;

const vector<int> DIRECTION_2 = {1, 3};

const vector<int> DIRECTION_3 = {1, 5, 17};

const vector<int> DIRECTION_4 = {1, 5, 23, 113};

const vector<vector<int>> K2DIRECTION = {{}, {}, DIRECTION_2, DIRECTION_3, DIRECTION_4};

constexpr int K2VERTEX_PRECOMP_FIRST[6] = {-1, -1, 4, 4, 4, 13};
constexpr int K2VERTEX_PRECOMP_LAST[6] = {-1, -1, 0, 1, 4, 4};

vector<int> cross_product_2(const vector<int>& x) {
    int c0 = x[1];
    int c1 = -x[0];

    return {c0, c1};
}

vector<int> cross_product_3(const vector<int>& x, const vector<int>& y) {
    int c0 = x[1] * y[2] - x[2] * y[1];

    int c1 = x[0] * y[2] - x[2] * y[0];
    c1 = -c1;

    int c2 = x[0] * y[1] - x[1] * y[0];

    return {c0, c1, c2};
}

vector<int> cross_product_4(const vector<int>& x, const vector<int>& y, const vector<int>& z) {
    int c0 = x[1] * y[2] * z[3] +
             y[1] * z[2] * x[3] +
             z[1] * x[2] * y[3] -
             z[1] * y[2] * x[3] -
             y[1] * x[2] * z[3] -
             x[1] * z[2] * y[3];

    int c1 = x[0] * y[2] * z[3] +
             y[0] * z[2] * x[3] +
             z[0] * x[2] * y[3] -
             z[0] * y[2] * x[3] -
             y[0] * x[2] * z[3] -
             x[0] * z[2] * y[3];
    c1 = -c1;

    int c2 = x[0] * y[1] * z[3] +
             y[0] * z[1] * x[3] +
             z[0] * x[1] * y[3] -
             z[0] * y[1] * x[3] -
             y[0] * x[1] * z[3] -
             x[0] * z[1] * y[3];

    int c3 = x[0] * y[1] * z[2] +
             y[0] * z[1] * x[2] +
             z[0] * x[1] * y[2] -
             z[0] * y[1] * x[2] -
             y[0] * x[1] * z[2] -
             x[0] * z[1] * y[2];

    c3 = -c3;
    return {c0, c1, c2, c3};
}

/*
 * Basically returns N choose K - 1 combinations.
 * The algorithm is taken from:
 * https://stackoverflow.com/questions/9430568/generating-combinations-in-c
 */
vector<vector<int>> generate_hyperplane_combinations(const vector<int>& hyperplanes, const int K) {
    int n = (int) hyperplanes.size();
    ASRTF(n >= K, "Should be a least k hyperplanes.");

    vector<bool> include_h(n);
    fill(include_h.begin(), include_h.begin() + K - 1, true);

    vector<vector<int>> combs;
    do {
        vector<int> comb(K - 1);
        int index = 0;
        for (int i = 0; i < n; i++) {
            if (include_h[i]) {
                comb[index] = hyperplanes[i];
                index++;
            }
        }
        // Can be performance critical, so normal assert.
        assert(index == K - 1 && "Expected that combination has size k - 1.");
        combs.push_back(comb);
    } while (prev_permutation(include_h.begin(), include_h.end()));
    ASRTF((int) combs.size() >= K, "Should be at least k combinations.");

    return combs;
}

struct Vertex {
    mpq_t* v{};
    // Contrary to the v, v_double doesn't have a leading 1, so it has K elements, not K + 1.
    vector<double> v_double;
    vector<int> incidence;
    int id{};
};

// TODO[gleb] Can also be used by comparing doubles before going to mpq - that can accelerate
// the comparison function.
struct MpqCompareVertex {
    bool operator() (const Vertex& lhs, const Vertex& rhs) const {
        // Can also compare with doubles and only then with
        const int K = (int) lhs.v_double.size();
        for (int i = 0; i < K; i++) {
            int cmp = mpq_cmp(lhs.v[i + 1], rhs.v[i + 1]);
            if (cmp < 0) {
                return true;
            } else if (cmp > 0) {
                return false;
            }
        }
        return false;
    }
};

inline int scalar_product(const vector<int>& coef1, const vector<int>& coef2) {
    // Can be performance critical, so normal assert.
    assert(coef1.size() == coef2.size() && "Vector sizes should match");

    int accum = 0;
    for (size_t i = 0; i < coef1.size(); i++) {
        accum += coef1[i] * coef2[i];
    }

    return accum;
}

bool vertices_divided_by_orthant(const mpq_t* v1, const mpq_t* v2, const int K) {
    for (int i = 0; i < K; i++) {
        int sign1 = mpq_sgn(v1[i + 1]);
        if (sign1 == 0) {
            continue;
        }
        int sign2 = mpq_sgn(v2[i + 1]);
        if (sign2 != 0 && sign1 != sign2) {
            return true;
        }
    }
    return false;
}

void print_Vertex(Vertex& vertex) {
    cout << "printing vertex" << endl;
    for (auto& v : vertex.v_double) {
        cout << v << " ";
    }
    cout << endl;
}

struct OctahedronToV_Helper {
    const int K;
    const int NUM_H;

    const vector<vector<int>>& COEFS;
    const vector<int>& MAIN_DIRECTION;
    const int PRECOMP_FIRST_SIZE;
    const int PRECOMP_LAST_SIZE;

    mpq_t* constraints {};
    vector<double> constraints_double;

    mpq_t* vertex_precomp_first {};
    mpq_t* vertex_precomp_last {};
    vector<double> vertex_precomp_first_double;
    vector<double> vertex_precomp_last_double;

    // Used to keep track of # of vertices processed so far. Is used to distribute ids to vertices.
    int vertex_count;

    set<Vertex, MpqCompareVertex> visited;
    set<Vertex, MpqCompareVertex> to_visit;
    vector<Adj> adjacencies;

    explicit OctahedronToV_Helper(const int K, const vector<double*>& A) :
        K(K),
        NUM_H(K2NUM_H[K]),
        COEFS(K2OCTAHEDRON_COEFS[K]),
        MAIN_DIRECTION(K2DIRECTION[K]),
        PRECOMP_FIRST_SIZE(K2VERTEX_PRECOMP_FIRST[K]),
        PRECOMP_LAST_SIZE(K2VERTEX_PRECOMP_LAST[K]),
        constraints_double(NUM_H),
        vertex_precomp_first_double(PRECOMP_FIRST_SIZE),
        vertex_precomp_last_double(PRECOMP_LAST_SIZE),
        vertex_count(0)
    {
        constraints = mpq_arr_create(NUM_H);
        for (int i = 0; i < NUM_H; i++) {
            mpq_set_d(constraints[i], A[i][0]);
            constraints_double[i] = A[i][0];
        }

        vertex_precomp_first = mpq_arr_create(PRECOMP_FIRST_SIZE);
        vertex_precomp_last = mpq_arr_create(PRECOMP_LAST_SIZE);
    }

    ~OctahedronToV_Helper() {
        // I don't free up dynamically allocated memory for vertices because these pointers are extracted and reused.
        mpq_arr_free(NUM_H, constraints);
        mpq_arr_free(PRECOMP_FIRST_SIZE, vertex_precomp_first);
        mpq_arr_free(PRECOMP_LAST_SIZE, vertex_precomp_last);
    }

    void do_precomputation_for_vertex(const Vertex& cur_vertex) {
        const mpq_t* cur_v = cur_vertex.v;
        const vector<double>& cur_v_double = cur_vertex.v_double;

        if (PRECOMP_FIRST_SIZE == 4) {
            vertex_precomp_first_double[0] = cur_v_double[0] + cur_v_double[1];
            vertex_precomp_first_double[1] = cur_v_double[0];
            vertex_precomp_first_double[2] = cur_v_double[0] - cur_v_double[1];
            vertex_precomp_first_double[3] = cur_v_double[1];

            mpq_add(vertex_precomp_first[0], cur_v[1], cur_v[2]);
            mpq_set(vertex_precomp_first[1], cur_v[1]);
            mpq_sub(vertex_precomp_first[2], cur_v[1], cur_v[2]);
            mpq_set(vertex_precomp_first[3], cur_v[2]);
        } else {
            vertex_precomp_first_double[0] = cur_v_double[0] + cur_v_double[1] + cur_v_double[2];
            vertex_precomp_first_double[1] = cur_v_double[0] + cur_v_double[1];
            vertex_precomp_first_double[2] = cur_v_double[0] + cur_v_double[1] - cur_v_double[2];
            vertex_precomp_first_double[3] = cur_v_double[0] + cur_v_double[2];
            vertex_precomp_first_double[4] = cur_v_double[0];
            vertex_precomp_first_double[5] = cur_v_double[0] - cur_v_double[2];
            vertex_precomp_first_double[6] = cur_v_double[0] - cur_v_double[1] + cur_v_double[2];
            vertex_precomp_first_double[7] = cur_v_double[0] - cur_v_double[1];
            vertex_precomp_first_double[8] = cur_v_double[0] - cur_v_double[1] - cur_v_double[2];
            vertex_precomp_first_double[9] = cur_v_double[1] + cur_v_double[2];
            vertex_precomp_first_double[10] = cur_v_double[1];
            vertex_precomp_first_double[11] = cur_v_double[1] - cur_v_double[2];
            vertex_precomp_first_double[12] = cur_v_double[2];

            // Slightly optimized version.
            mpq_add(vertex_precomp_first[1], cur_v[1], cur_v[2]);
            mpq_add(vertex_precomp_first[0], vertex_precomp_first[1], cur_v[3]);
            mpq_sub(vertex_precomp_first[2], vertex_precomp_first[1], cur_v[3]);
            mpq_add(vertex_precomp_first[3], cur_v[1], cur_v[3]);
            mpq_set(vertex_precomp_first[4], cur_v[1]);
            mpq_sub(vertex_precomp_first[5], cur_v[1], cur_v[3]);
            mpq_sub(vertex_precomp_first[7], cur_v[1], cur_v[2]);
            mpq_add(vertex_precomp_first[6], vertex_precomp_first[7], cur_v[3]);
            mpq_sub(vertex_precomp_first[8], vertex_precomp_first[7], cur_v[3]);
            mpq_add(vertex_precomp_first[9], cur_v[2], cur_v[3]);
            mpq_set(vertex_precomp_first[10], cur_v[2]);
            mpq_sub(vertex_precomp_first[11], cur_v[2], cur_v[3]);
            mpq_set(vertex_precomp_first[12], cur_v[3]);

        }

        if (PRECOMP_LAST_SIZE == 1) {
            vertex_precomp_last_double[0] = cur_v_double[K-1];

            mpq_set(vertex_precomp_last[0], cur_v[K]);
        } else if (PRECOMP_LAST_SIZE == 4) {
            vertex_precomp_last_double[0] = cur_v_double[K-2] + cur_v_double[K-1];
            vertex_precomp_last_double[1] = cur_v_double[K-2];
            vertex_precomp_last_double[2] = cur_v_double[K-2] - cur_v_double[K-1];
            vertex_precomp_last_double[3] = cur_v_double[K-1];

            mpq_add(vertex_precomp_last[0], cur_v[K-1], cur_v[K]);
            mpq_set(vertex_precomp_last[1], cur_v[K-1]);
            mpq_sub(vertex_precomp_last[2], cur_v[K-1], cur_v[K]);
            mpq_set(vertex_precomp_last[3], cur_v[K]);
        } else {
            ASRTF(PRECOMP_LAST_SIZE == 0, "Expected size 0.");
        }
    }

    vector<vector<int>> generate_directions(const Vertex& vertex) {
        vector<vector<int>> hyperplane_combs = generate_hyperplane_combinations(vertex.incidence, K);

        set<vector<int>> directions;
        for (const auto& comb : hyperplane_combs) {
            vector<int> direction;
            if (K == 2) {
                direction = cross_product_2(COEFS[comb[0]]);
            } else if (K == 3) {
                direction = cross_product_3(COEFS[comb[0]], COEFS[comb[1]]);
            } else {
                direction = cross_product_4(COEFS[comb[0]], COEFS[comb[1]], COEFS[comb[2]]);
            }

            // It is possible that a combination of hyperplanes is linearly dependent and doesn't create a
            // direction - in that case I will just skip it.
            bool at_least_one_not_zero = false;
            for (int i = 0; i < K; i++) {
                if (direction[i] != 0) {
                    at_least_one_not_zero = true;
                    break;
                }
            }
            if (at_least_one_not_zero) {
                directions.insert(direction);
            }
        }

        // Now for every direction we determine whether it's positive, negative or non-viable.
        // There are 4 possible cases:
        // (1) The direction is positive iff scalar products with all hyperplanes
        //      give >= 0 and gives > 0 with at least 1 hyperplane.
        // (2) The direction is negative iff -//- but with <= 0 and < 0.
        // (3) The direction is non-viable iff there are hyperplanes with which scalar
        //      product gives < 0 and hyperplanes with which scalar product gives > 0.
        // (4) There is an error if for some direction scalar product with all hyperplanes is 0.
        //      It means that the vertex was not extreme in the first place.

        // It is important that it is a set - otherwise there is a chance that
        // direction is negative, I will inverse it and will create a duplicate.
        set<vector<int>> selected_directions;

        // Intentionally making a copy so I can modify it.
        for (auto dir : directions) {
            bool is_positive = false;
            bool is_negative = false;

            for (int hi : vertex.incidence) {
                const vector<int>& h = COEFS[hi];

                int product = scalar_product(dir, h);
                if (product > 0) {
                    is_positive = true;
                    if (is_negative) {
                        // when is_negative and is_positive the direction is non-viable, I can stop early.
                        break;
                    }
                } else if (product < 0) {
                    is_negative = true;
                    if (is_positive) {
                        // when is_negative and is_positive the direction is non-viable, I can stop early.
                        break;
                    }
                }
            }

            ASRTF(is_positive || is_negative,
                  "Going in any direction doesn't violate any hyperplanes. The vertex is not extreme.");
            if (is_positive && is_negative) {
                // The direction is non-viable - going in any direction will violate at least one of the incident
                // hyperplanes.
                continue;
            }
            if (is_negative) {
                // To keep code simple, I always want to go in positive direction relative to the main direction.
                for (int i = 0; i < K; i++) {
                    dir[i] = -dir[i];
                }
            }
            int product_with_main_direction = scalar_product(MAIN_DIRECTION, dir);
            // The product with the main direction can't be zero, because I select main direction in such way that no
            // edge of the octahedra can be parallel to this edge. To prove it, I try all various combinations of
            // hyperplanes to produce the edge and make sure that the scalar product is not zero.
            ASRTF(product_with_main_direction != 0, "Scalar product with the direction can't be zero.");
            if (product_with_main_direction > 0) {
                selected_directions.insert(dir);
            }
        }

        return vector<vector<int>>(selected_directions.begin(), selected_directions.end());
    }

    void visit_vertex() {
        ASRTF(!to_visit.empty(), "No vertices left to visit.");
//        cout << "visiting vertex" << endl;
        auto vertex_iter = to_visit.begin();
        Vertex vertex = *vertex_iter;

        to_visit.erase(vertex_iter);
        visited.insert(vertex);

        vector<vector<int>> directions = generate_directions(vertex);
//        cout << "directions size " << directions.size() << endl;

//        for (auto& dir : directions) {
//            cout << "print direction" << endl;
//            for (auto& el : dir) {
//                cout << el << " ";
//            }
//            cout << endl;
//        }

        do_precomputation_for_vertex(vertex);

        vector<int> denom_all(directions.size());
        vector<double> closest_t_double_all(directions.size());
        mpq_t* closest_t_all = mpq_arr_create(directions.size());
        // The hyperplanes that will be incident to the vertex reached by new direction.
        vector<vector<int>> new_hyperplanes_all(directions.size());

        mpq_t cur_t, numer, denom;
        mpq_inits(cur_t, numer, denom, NULL);

        set_t incident = set_create(NUM_H);
        for (int hi : vertex.incidence) {
            set_enable_bit(incident, hi);
        }
        for (int hi = 0; hi < NUM_H; hi++) {
            if (set_test_bit(incident, hi)) {
                // If the vertex is incident to the current vertex then obviously t = 0.
                continue;
            }
            set_t skip_direction = set_create(directions.size());
            const vector<int>& h = COEFS[hi];

            int precomp_first_index = -1;
            int precomp_last_index = -1;

            if (PRECOMP_FIRST_SIZE == 4) {
                precomp_first_index = 3 * (1 - h[0]) + 1 - h[1];
            } else {
                precomp_first_index = 9 * (1 - h[0]) + 3 * (1 - h[1]) + 1 - h[2];
            }

            if (PRECOMP_LAST_SIZE == 1) {
                precomp_last_index = 1 - h[K - 1];
            } else if (PRECOMP_LAST_SIZE == 4) {
                precomp_last_index = 3 * (1 - h[K - 2]) + 1 - h[K - 1];
            }

            double numer_double = constraints_double[hi];

            if (precomp_first_index < PRECOMP_FIRST_SIZE) {
                numer_double += vertex_precomp_first_double[precomp_first_index];
            } else if (precomp_first_index > PRECOMP_FIRST_SIZE) {
                numer_double -= vertex_precomp_first_double[2 * PRECOMP_FIRST_SIZE - precomp_first_index];
            }

            if (PRECOMP_LAST_SIZE != 0) {
                if (precomp_last_index < PRECOMP_LAST_SIZE) {
                    numer_double += vertex_precomp_last_double[precomp_last_index];
                } else if (precomp_last_index > PRECOMP_LAST_SIZE) {
                    numer_double -= vertex_precomp_last_double[2 * PRECOMP_LAST_SIZE - precomp_last_index];
                }
            }

            bool at_least_one_left_to_process = false;
            for (int dir_i = 0; dir_i < (int) directions.size(); dir_i++) {
                // I negate denominator because I have numerator in negated form.
                int denom_int = -scalar_product(directions[dir_i], h);
                if (denom_int == 0) {
                    // The hyperplane is parallel to this direction.
                    set_enable_bit(skip_direction, dir_i);
                    continue;
                }
                double cur_t_double = numer_double / (double) denom_int;
                if (cur_t_double < -DOUBLE_PRECISION) {
                    // Rejecting the hyperplane early based on the floating point
                    // arithmetic because hyperplane intersects in the negative direction.
                    set_enable_bit(skip_direction, dir_i);
                    continue;
                }
                if (!new_hyperplanes_all[dir_i].empty() &&
                    cur_t_double > closest_t_double_all[dir_i] + DOUBLE_PRECISION) {
                    // Rejecting hyperplane because cur_t_double is significantly
                    // bigger than best result.
                    set_enable_bit(skip_direction, dir_i);
                    continue;
                }
                at_least_one_left_to_process = true;
                denom_all[dir_i] = denom_int;
            }

            if (!at_least_one_left_to_process) {
                // Based on floating point arithmetics this hyperplane can't give next
                // vertex for any of the directions.
                // Thus just moving forward to the next hyperplane.
                continue;
            }

            mpq_set(numer, constraints[hi]);

            if (precomp_first_index < PRECOMP_FIRST_SIZE) {
                mpq_add(numer, numer, vertex_precomp_first[precomp_first_index]);
            } else if (precomp_first_index > PRECOMP_FIRST_SIZE) {
                mpq_sub(numer, numer, vertex_precomp_first[2 * PRECOMP_FIRST_SIZE - precomp_first_index]);
            }

            if (PRECOMP_LAST_SIZE != 0) {
                if (precomp_last_index < PRECOMP_LAST_SIZE) {
                    mpq_add(numer, numer, vertex_precomp_last[precomp_last_index]);
                } else if (precomp_last_index > PRECOMP_LAST_SIZE) {
                    mpq_sub(numer, numer, vertex_precomp_last[2 * PRECOMP_LAST_SIZE - precomp_last_index]);
                }
            }

            int sign = mpq_sgn(numer);
            // Can be performance critical, so normal assert.
            assert(sign != 0 &&
                   "Sign can't be zero, because the hyperplanes incident to the current vertex were rejected earlier.");

            for (int dir_i = 0; dir_i < (int) directions.size(); dir_i++) {
                if (set_test_bit(skip_direction, dir_i)) {
                    // This direction was rejected earlier based on floating point arithmetic.
                    continue;
                }

                int denom_int = denom_all[dir_i];
                if ((denom_int < 0 && sign > 0) || (denom_int > 0 && sign < 0)) {
                    // cur_t will be < 0 - this means that this direction is of no interest
                    // because I'm only looking at the positive directions.
                    continue;
                }

                mpq_set_si(denom, denom_int, 1);
                mpq_div(cur_t, numer, denom);

                vector<int>& new_hyperplanes = new_hyperplanes_all[dir_i];
                if (new_hyperplanes.empty()) {
//                    cout << "setting new cur_t " << mpq_get_d(cur_t) << endl;
                    new_hyperplanes.push_back(hi);
                    mpq_set(closest_t_all[dir_i], cur_t);
                    closest_t_double_all[dir_i] = mpq_get_d(closest_t_all[dir_i]);
                    continue;
                }
                int comp = mpq_cmp(cur_t, closest_t_all[dir_i]);
                if (comp == 0) {
                    new_hyperplanes.push_back(hi);
                } else if (comp < 0) {
                    new_hyperplanes = {hi};
                    mpq_set(closest_t_all[dir_i], cur_t);
                    closest_t_double_all[dir_i] = mpq_get_d(closest_t_all[dir_i]);
                }
            }
            set_free(skip_direction);
        }

        set_free(incident);
        mpq_clears(cur_t, numer, denom, NULL);

        mpq_t dir_elem;
        mpq_init(dir_elem);

        for (size_t dir_i = 0; dir_i < directions.size(); dir_i++) {
            const vector<int>& dir = directions[dir_i];
            vector<int>& new_hyperplanes = new_hyperplanes_all[dir_i];
            ASRTF(!new_hyperplanes.empty(),
                  "At least one hyperplane should intersect direction. Otherwise octahedron is unbounded.");

            mpq_t* new_v = mpq_arr_create(K + 1);
            mpq_set_si(new_v[0], 1, 1);

            vector<double> new_v_double(K);
            mpq_t& closest_t = closest_t_all[dir_i];
            for (int i = 0; i < K; i++) {
                if (dir[i] == 0) {
                    mpq_set(new_v[i + 1], vertex.v[i + 1]);
                    new_v_double[i] = vertex.v_double[i];
                    continue;
                }
                mpq_set_si(dir_elem, dir[i], 1);
                mpq_mul(dir_elem, dir_elem, closest_t);
                mpq_add(new_v[i + 1], vertex.v[i + 1], dir_elem);
                new_v_double[i] = mpq_get_d(new_v[i + 1]);
            }

            // Optimistically create new_vertex. It might have already been computed.
            // In that case we should just skip this vertex.
            Vertex new_vertex {new_v, new_v_double, new_hyperplanes, vertex_count};
            bool save_adjacency = vertices_divided_by_orthant(vertex.v, new_v, K);
//            cout << "here save adjacency is " << save_adjacency << endl;

            // TODO[gleb] I can optimize this code by using incidence as a vertex identifier.
            // In that case I will not need to create final vertex and very quickly check for existence.
            auto search_visited = visited.find(new_vertex);
            if (search_visited != visited.end()) {
                if (save_adjacency) {
                    adjacencies.emplace_back(min(search_visited->id, vertex.id), max(search_visited->id, vertex.id));
                }
                mpq_arr_free(K + 1, new_v);
                continue;
            }

            auto search_to_visit = to_visit.find(new_vertex);
            if (search_to_visit != to_visit.end()) {
                if (save_adjacency) {
                    adjacencies.emplace_back(min(search_to_visit->id, vertex.id), max(search_to_visit->id, vertex.id));
                }
                mpq_arr_free(K + 1, new_v);
                continue;
            }

            // Will actually be adding this vertex, thus increase the count.
            if (save_adjacency) {
                adjacencies.emplace_back(min(new_vertex.id, vertex.id), max(new_vertex.id, vertex.id));
            }

            // Hyperplanes that were incident to the current vertex and
            // parallel to the direction will also be incident to
            // new vertex.
            for (int hi : vertex.incidence) {
                int product = scalar_product(dir, COEFS[hi]);

                if (product == 0) {
                    new_vertex.incidence.push_back(hi);
                }
            }

            // Have added the vertex, thus increasing the count.
            vertex_count++;
            to_visit.insert(new_vertex);
        }

        mpq_clear(dir_elem);
        mpq_arr_free(directions.size(), closest_t_all);
    }

    Vertex get_first_vertex() {
        dd_MatrixPtr lp_mat = dd_CreateMatrix(NUM_H, K + 1);
        for (int hi = 0; hi < NUM_H; hi++) {
            mpq_t* row = lp_mat->matrix[hi];
            const vector<int>& coef = COEFS[hi];
            mpq_set(row[0], constraints[hi]);
            for (int j = 0; j < K; j++) {
                mpq_set_si(row[j + 1], coef[j], 1);
            }
        }

        mpq_set_si(lp_mat->rowvec[0], 0, 1);
        for (int i = 0; i < K; i++) {
            mpq_set_si(lp_mat->rowvec[i + 1], MAIN_DIRECTION[i], 1);
        }

        dd_ErrorType error = dd_NoError;

        dd_LPPtr lp = dd_Matrix2LP(lp_mat, &error);
        lp->objective = dd_LPmin;
        ASRTF(error == dd_NoError, "Creating LP failed.");

        dd_LPSolve(lp, dd_DualSimplex, &error);
        ASRTF(error == dd_NoError, "Solving LP failed.");

        mpq_t* found_vertex = lp->sol;
        ASRTF(!mpq_cmp_si(found_vertex[0], 1, 1), "First element of the found vertex should be 1.");

        mpq_t* first_vertex = mpq_arr_create(K + 1);
        vector<double> first_vertex_double(K);
        mpq_set_si(first_vertex[0], 1, 1);
        for (int i = 0; i < K; i++) {
            mpq_set(first_vertex[i + 1], found_vertex[i + 1]);
            first_vertex_double[i] = mpq_get_d(found_vertex[i + 1]);
        }

        dd_FreeMatrix(lp_mat);
        dd_FreeLPData(lp);

        vector<int> first_inc;

        mpq_t accumulator;
        mpq_init(accumulator);
        for (int hi = 0; hi < NUM_H; hi++) {
            mpq_set_si(accumulator, 0, 1);
            const vector<int>& h_coef = COEFS[hi];
            for (int j = 0; j < K; j++) {
                if (h_coef[j] == 1) {
                    mpq_sub(accumulator, accumulator, first_vertex[j + 1]);
                } else if (h_coef[j] == -1) {
                    mpq_add(accumulator, accumulator, first_vertex[j + 1]);
                }
            }
            int comparison = mpq_cmp(constraints[hi], accumulator);
            ASRTF(comparison >= 0, "First vertex violates some constraint.");
            if (comparison == 0) {
                first_inc.push_back(hi);
            }
        }
        mpq_clear(accumulator);

        ASRTF((int) first_inc.size() >= K, "The vertex should be incident to at least k constraints.");

        vertex_count++;
        return {first_vertex, first_vertex_double, first_inc, 0};
    }

    void compute() {
        Vertex first_vertex = get_first_vertex();

        to_visit.insert(first_vertex);

        while (!to_visit.empty()) {
            visit_vertex();
        }

        ASRTF(vertex_count == (int) visited.size(), "Vertex id should equal the number of visited vertexes.");
    }

    OctahedronV extract_OctahedronV() {
        size_t num_vertices = visited.size();
        vector<mpq_t*> V(num_vertices);
        vector<set_t> incidence = set_arr_create(num_vertices, NUM_H);

        for (auto& vertex : visited) {
            int id = vertex.id;
            ASRTF(0 <= id && id < (int) num_vertices, "Id is out of range.");
            V[id] = vertex.v;
            auto& inc = incidence[id];

            for (int hi : vertex.incidence) {
                set_enable_bit(inc, hi);
            }
        }

        return {K + 1, V, adjacencies, incidence};
    }
};

OctahedronV compute_octahedron_V(const int K, const vector<double*>& A) {
    ASRTF(2 <= K && K <= 4, "K should be within allowed range.");
    ASRTF((int) A.size() == K2NUM_H[K], "Unexpected number of rows in the input.");
    OctahedronToV_Helper computation(K, A);
    computation.compute();
    return computation.extract_OctahedronV();
}

OctahedronV compute_V_with_cdd(const int K, const vector<double*>& A) {
    ASRTF(2 <= K, "K should be within allowed range.");

    dd_MatrixPtr cdd_hrep = fp_mat_to_cdd(K + 1, A);
    cdd_hrep->representation = dd_Inequality;
    dd_ErrorType err = dd_NoError;

    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(cdd_hrep, &err);
    ASRTF(err == dd_NoError, "Converting matrix to polytope failed with error " + to_string(err));
    dd_MatrixPtr cdd_vrep = dd_CopyGenerators(poly);

    int num_v = cdd_vrep->rowsize;

    ASRTF(num_v > 0, "V-representation is expected to have at least 1 vertex.");
    for (int i = 0; i < num_v; i++) {
        ASRTF(!mpq_cmp_si(cdd_vrep->matrix[i][0], 1, 1), "Polytope is not bounded.");
    }

    vector<mpq_t*> V(num_v);
    for (int i = 0; i < num_v; i++) {
        V[i] = mpq_arr_copy(K + 1, cdd_vrep->matrix[i]);
    }

    dd_SetFamilyPtr cdd_adjacencies = dd_CopyAdjacency(poly);
    dd_SetFamilyPtr cdd_incidence = dd_CopyIncidence(poly);

    vector<Adj> adjacencies;
    for (int vi = 0; vi < num_v; vi++) {
        const mpq_t* v = V[vi];
        vector<int> signs(K + 1);
        for (int d = 1; d <= K; d++) {
            signs[d] = mpq_sgn(v[d]);
        }
        set_type set_i = cdd_adjacencies->set[vi];
        for (int vj = vi + 1; vj < num_v; vj++) {
            // CDD's sets are 1-based.
            if (!set_member(vj + 1, set_i)) {
                // Checking if vi and vj are adjacent at all.
                continue;
            }
            if (vertices_divided_by_orthant(v, V[vj], K)) {
                adjacencies.emplace_back(vi, vj);
            }
        }
    }

    vector<set_t> incidence(num_v);
    for (int v = 0; v < num_v; v++) {
        incidence[v] = set_from_cdd(cdd_incidence->set[v]);
    }

    dd_FreeMatrix(cdd_hrep);
    dd_FreeMatrix(cdd_vrep);
    dd_FreePolyhedra(poly);
    dd_FreeSetFamily(cdd_adjacencies);
    dd_FreeSetFamily(cdd_incidence);

    return {K + 1, V, adjacencies, incidence};
}

map<Quadrant, QuadrantInfo> compute_quadrants_with_cdd(const int K, const vector<double*>& A) {
    ASRTF(1 <= K && K <= 5, "K should be within allowed range.");

    vector<Quadrant> quadrants {{}};
    for (int i = 0; i < K; i++) {
        vector<Quadrant> next_quadrants;
        next_quadrants.reserve(quadrants.size() * 2);

        for (auto& quadrant : quadrants) {
            quadrant.push_back(MINUS);
            next_quadrants.push_back(quadrant);
            quadrant.back() = PLUS;
            next_quadrants.push_back(quadrant);
        }
        quadrants = next_quadrants;
    }
    ASRTF((int) quadrants.size() == POW2[K], "Sanity checking number of quadrants.");

    const int NUM_H = (int) A.size();

    dd_MatrixPtr cdd_A = dd_CreateMatrix(NUM_H + K, K + 1);
    cdd_A->representation = dd_Inequality;

    for (int i = 0; i < NUM_H; i++) {
        mpq_t* cdd_row = cdd_A->matrix[i];
        const double* row = A[i];
        for (int j = 0; j < K + 1; j++) {
            mpq_set_d(cdd_row[j], row[j]);
        }
    }

    map<Quadrant, QuadrantInfo> quadrant2info;
    for (const auto& quadrant : quadrants) {
        for (int xi = 0; xi < K; xi++) {
            mpq_t* row = cdd_A->matrix[NUM_H + xi];
            mpq_arr_set_zero(K + 1, row);
            if (quadrant[xi] == MINUS) {
                mpq_set_si(row[xi + 1], -1, 1);
            } else {
                mpq_set_si(row[xi + 1], 1, 1);
            }
        }
        dd_ErrorType err = dd_NoError;
        dd_PolyhedraPtr poly = dd_DDMatrix2Poly(cdd_A, &err);

        ASRTF(err == dd_NoError, "Converting matrix to polytope failed with error " + to_string(err));

        dd_MatrixPtr cdd_V = dd_CopyGenerators(poly);
        dd_SetFamilyPtr cdd_incidence = dd_CopyIncidence(poly);
        const size_t num_v = cdd_V->rowsize;
        ASRTF(
                cdd_incidence->famsize == (int) num_v && \
                cdd_incidence->setsize == NUM_H + K + 1,
                "Sanity checking cdd_incidence size.");

        vector<mpq_t*> V(num_v);
        // Here H include the orthant.
        vector<set_t> incidence(num_v);

        for (size_t v = 0; v < num_v; v++) {
            ASRTF(!mpq_cmp_si(cdd_V->matrix[v][0], 1, 1), "Polytope is not bounded.");
            V[v] = mpq_arr_copy(K + 1, cdd_V->matrix[v]);
            incidence[v] = set_from_cdd(cdd_incidence->set[v]);
        }
        dd_FreePolyhedra(poly);
        dd_FreeMatrix(cdd_V);
        dd_FreeSetFamily(cdd_incidence);

        quadrant2info[quadrant] = {K + 1, V, incidence};
    }

    dd_FreeMatrix(cdd_A);
    return quadrant2info;
}
