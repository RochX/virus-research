#include "std_vector_functions.hpp"

// need to explicitly declare what types are used with the template
template void std_vector_functions::append_vector(std::vector<EigenType::Vector6f>& v1, std::vector<EigenType::Vector6f>& v2, bool removeDups);
template void std_vector_functions::push_backIfNotInVector(std::vector<float> &vector, float element, float epsilon);

// append the contents of v2 to v1, possibly checking for duplicates
template <typename T>
void std_vector_functions::append_vector(std::vector<T>& v1, std::vector<T>& v2, bool removeDups) {
    for (const T &v : v2) {
        // if bool is true only append if it is not a duplicate
        if (removeDups) {
            if (std::find(v1.begin(), v1.end(), v) == v1.end()) {
                v1.push_back(v);
            }
        }
        // otherwise just append without checking
        else {
            v1.push_back(v);
        }
    }
}

template <typename T>
void std_vector_functions::push_backIfNotInVector(std::vector<T> &vector, T element, T epsilon) {
    bool addElement = true;
    for (T t : vector) {
        if (std::fabs(t-element) <= epsilon) {
            addElement = false;
            break;
        }
    }

    if (addElement)
        vector.push_back(element);
}