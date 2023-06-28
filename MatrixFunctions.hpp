#include <iostream>
#include <vector>
#include <Eigen/Core>

#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_MATRIXFUNCTIONS_HPP
#define VIRUS_RESEARCH_MATRIXFUNCTIONS_HPP

namespace MatrixFunctions {
    void print_size(const Eigen::MatrixXf& b);
    std::vector<float> entriesOfMatrix(const Eigen::MatrixXf& m, bool remove_duplicates = false);
    bool entriesOfMatrixAreOfParticularValues(const Eigen::MatrixXf& matrix, const std::vector<float>& values);
    void fixZeroEntries(Eigen::MatrixXf& matrix);
    void fixZeroEntries(EigenType::Matrix6f& matrix);
    void fixZeroEntries(EigenType::Matrix6fx3f& matrix);
    bool matrixIsScalarOfIdentity(const EigenType::Matrix6f& matrix);
    bool matrixIsPlusMinusIdentity(const EigenType::Matrix6f& matrix);
} // MatrixFunctions

#endif //VIRUS_RESEARCH_MATRIXFUNCTIONS_HPP
