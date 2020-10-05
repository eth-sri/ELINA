#pragma once

#include <string>
#include <stdexcept>

using namespace std;

#define ASRTF(A,B) asrt(A,B, string(__FILE__) + " " + __FUNCTION__ + " " + to_string(__LINE__))

// Inlining this function makes a big performance difference.
inline void asrt(bool condition, const string &error_message = "", const string &func_name = "") {
    if (!condition) {
        if (func_name.empty()) {
            throw runtime_error(error_message);
        } else {
            throw runtime_error(func_name + " : " + error_message);
        }
    }
}
