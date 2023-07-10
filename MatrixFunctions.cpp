#include "MatrixFunctions.hpp"

namespace MatrixFunctions {
    void print_size(const Eigen::MatrixXf& b)
    {
        std::cout << "size (rows, cols): " << b.size() << " (" << b.rows()
                  << ", " << b.cols() << ")" << std::endl;
    }

    std::vector<float> entriesOfMatrix(const Eigen::MatrixXf& m, bool remove_duplicates) {
        std::vector<float> entries;
        float curr;
        for (size_t i = 0; i < m.size(); i++) {
            curr = *(m.data()+i);
            if (remove_duplicates) {
                bool found = [&] {
                    for (float f : entries) {
                        if (float_functions::FloatsAreApproxEqual(curr, f))
                            return true;
                    }
                    return false;
                }();
                if (!found)
                    entries.push_back(curr);
            }
            else
                entries.push_back(curr);
        }
        return entries;
    }

    bool entriesOfMatrixAreOfParticularValues(const Eigen::MatrixXf& matrix, const std::vector<float>& values) {
        bool entry_is_valid;
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                entry_is_valid = false;

                // check if it is any of the valid values
                for (float val : values) {
                    if (float_functions::FloatsAreApproxEqual(matrix(i,j), val)) {
                        entry_is_valid = true;
                        break;
                    }
                }

                if (!entry_is_valid)
                    return false;
            }
        }

        // every entry is valid
        return true;
    }

    void fixZeroEntries(Eigen::MatrixXf& matrix) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                if (float_functions::FloatsAreApproxEqual(matrix(i,j), 0)) {
                    matrix(i,j) = 0;
                }
            }
        }
    }

    void fixZeroEntries(EigenType::Matrix6f& matrix) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                if (float_functions::FloatIsApproxZero(matrix(i,j))) {
                    matrix(i,j) = 0;
                }
            }
        }
    }

    void fixZeroEntries(EigenType::Matrix6fx3f& matrix) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                if (float_functions::FloatIsApproxZero(matrix(i,j))) {
                    matrix(i,j) = 0;
                }
            }
        }
    }

    bool matrixIsScalarOfIdentity(const EigenType::Matrix6f& matrix) {
        return (matrix * 1.0/matrix(0,0)).isApprox(EigenType::Matrix6f::Identity());
    }

    bool matrixIsPlusMinusIdentity(const EigenType::Matrix6f& matrix) {
        return matrix.isApprox(EigenType::Matrix6f::Identity()) || (-1*matrix).isApprox(EigenType::Matrix6f::Identity());
    }
} // MatrixFunctions