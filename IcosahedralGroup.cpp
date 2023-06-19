#include "IcosahedralGroup.hpp"

// constructor generates the group
IcosahedralGroup::IcosahedralGroup() {
    // A has order 2
    generatorA << -1, 0, 0, 0, 0, 0,
            0, -1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0;

    // B has order 3
    generatorB << 0, -1, 0, 0, 0, 0,
            0, 0, -1, 0, 0, 0,
            1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0;

    createGroupFromTwoGeneratorsAndKnownSize(generatorA, generatorB, 60);
}

EigenType::Matrix6f IcosahedralGroup::matrixFormOfCentralizer(float z, float x) {
    EigenType::Matrix6f c;
    c << z,x,-x,-x,x,x,
            x,z,x,-x,-x,x,
            -x,x,z,x,-x,x,
            -x,-x,x,z,x,x,
            x,-x,-x,x,z,x,
            x,x,x,x,x,z;
    return c;
}

bool IcosahedralGroup::checkIfInCentralizer(EigenType::Matrix6f m) {
    return m.isApprox(matrixFormOfCentralizer(m(0,0), m(0,1)));
}