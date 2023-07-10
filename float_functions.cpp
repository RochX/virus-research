#include "float_functions.hpp"

namespace float_functions {
    bool FloatIsApproxZero(float x) {
        const float epsilon = 0.00001;
        return std::fabs(x) < epsilon;
    }

    bool FloatsAreApproxEqual(float x, float y) {
        if (std::fabs(x) < 0.00001)
            return FloatIsApproxZero(y);
        if (std::fabs(y) < 0.00001)
            return FloatIsApproxZero(x);

        const float relative_difference_factor = 0.0001;    // 0.01%
        const float greater_magnitude = std::max(std::fabs(x),std::fabs(y));

        return std::fabs(x-y) < relative_difference_factor * greater_magnitude;
    }
} // float_functions