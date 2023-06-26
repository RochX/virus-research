#include "std_vector_functions.hpp"

// need to explicitly declare what types are used with the template
template void std_vector_functions::append_vector(std::vector<float>& v1, std::vector<float>& v2, bool removeDups);
template void std_vector_functions::append_vector(std::vector<EigenType::Vector6f>& v1, std::vector<EigenType::Vector6f>& v2, bool removeDups);
template void std_vector_functions::append_vector(std::vector<EigenType::Matrix6f>& v1, std::vector<EigenType::Matrix6f>& v2, bool removeDups);
template void std_vector_functions::push_backIfNotInVector(std::vector<float> &vector, float element, float epsilon);
template std::vector<EigenType::Vector6f> std_vector_functions::unravelTwoDimVector(std::vector<std::vector<EigenType::Vector6f>> two_dim_vector, bool remove_duplicates);

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

template <typename T>
std::vector<T> std_vector_functions::unravelTwoDimVector(std::vector<std::vector<T>> two_dim_vector, bool remove_duplicates) {
    std::vector<T> unraveled;

    for (std::vector<T> vector : two_dim_vector) {
        for (T element : vector) {
            if (remove_duplicates) {
                // if element not in unraveled version, add it
                if (std::find(unraveled.begin(), unraveled.end(), element) == unraveled.end()) {
                    unraveled.push_back(element);
                }
            }
            // just push back all elements if not checking for duplicates
            else {
                unraveled.push_back(element);
            }
        }
    }

    return unraveled;
}