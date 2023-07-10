#include <vector>
#include <cmath>
#include <algorithm>

#include "EigenTypes.hpp"
#include "float_functions.hpp"

#ifndef VIRUS_RESEARCH_STD_VECTOR_FUNCTIONS_HPP
#define VIRUS_RESEARCH_STD_VECTOR_FUNCTIONS_HPP

namespace std_vector_functions {
    template <typename T>
    void append_vector(std::vector<T>& v1, std::vector<T>& v2, bool removeDups);
    void append_vector(std::vector<float> &v1, std::vector<float> &v2, bool removeDups);

    template <typename T>
    void push_backIfNotInVector(std::vector<T> &vector, T element, T epsilon = 0);

    template <typename T>
    std::vector<T> unravelTwoDimVector(std::vector<std::vector<T>> two_dim_vector, bool remove_duplicates = false);
}

#endif //VIRUS_RESEARCH_STD_VECTOR_FUNCTIONS_HPP
