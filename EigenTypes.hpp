//
// Created by xavier on 6/16/2023.
//
#include <Eigen/Dense>

#ifndef VIRUS_RESEARCH_MATRIXTYPE_H
#define VIRUS_RESEARCH_MATRIXTYPE_H

namespace EigenType {
    typedef Eigen::Vector<float, 6> Vector6f;
    typedef Eigen::Matrix<float, 6, 6> Matrix6f;
    typedef Eigen::Matrix<float, 6, 3> Matrix6fx3f;
    static Eigen::IOFormat COMMA_SEP_VALS(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
    static Eigen::IOFormat TAB_INDENT(Eigen::StreamPrecision, 0, " ", "\n", "\t");
}

#endif //VIRUS_RESEARCH_MATRIXTYPE_H
