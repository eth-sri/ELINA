#pragma once

#include <string>
#include <stdexcept>

#define ASRTF(A,B) asrt(A,B, std::string(__FILE__) + " " + __FUNCTION__ + " " + std::to_string(__LINE__))

// Inlining this function makes a big performance difference.
inline void asrt(bool condition, const std::string &error_message = "", const std::string &func_name = "") {
    if (!condition) {
        if (func_name.empty()) {
            throw std::runtime_error(error_message);
        } else {
            throw std::runtime_error(func_name + " : " + error_message);
        }
    }
}
