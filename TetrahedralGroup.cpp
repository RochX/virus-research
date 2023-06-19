//
// Created by xavier on 6/19/2023.
//

#include "TetrahedralGroup.hpp"

TetrahedralGroup::TetrahedralGroup() {
    // A has order 2
    generatorA << 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, -1, 0, 0,
        0, 0, 0, 0, -1, 0,
        0, 1, 0, 0, 0, 0;

    // B has order 3
    generatorB << 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0,
        0, -1, 0, 0, 0, 0,
        0, 0, -1, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0;

    createGroupFromTwoGeneratorsAndKnownSize(generatorA, generatorB, 12);
}

EigenType::Matrix6f TetrahedralGroup::matrixFormOfCentralizer(float x, float y, float z, float t) {
    EigenType::Matrix6f c;

    c << z, -x, -y, -t, t, -x,
        t, z, t, x, x, y,
        -y, -x, z, t, -t, -x,
        x, -t, -x, z, y, t,
        -x, -t, x, y, z, t,
        t, y, t, -x, -x, z;

    return c;

}

bool TetrahedralGroup::checkIfInCentralizer(EigenType::Matrix6f m) {
    return m.isApprox(matrixFormOfCentralizer(m(4, 2), m(4, 3), m(4, 4), m(4, 5)));
}


