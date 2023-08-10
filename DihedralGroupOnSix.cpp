#include "DihedralGroupOnSix.hpp"

DihedralGroupOnSix::DihedralGroupOnSix() {
    // A has order 3
    generatorA << 0, 0, 0, 0, 0, 1,
                0, 0, 0, 1, 0, 0,
                0, -1, 0, 0, 0, 0,
                0, 0, -1, 0, 0, 0,
                1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0;

    // B has order 2
    generatorB << 0, 0, 0, 0, -1, 0,
                0, 0, 0, -1, 0, 0,
                0, 0, -1, 0, 0, 0,
                0, -1, 0, 0, 0, 0,
                -1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, -1;

    // D_10 has 10 elements
    createGroupFromTwoGeneratorsAndKnownSize(generatorA, generatorB, 6);
}

EigenType::Matrix6f DihedralGroupOnSix::matrixFormOfCentralizer(float x, float y, float z, float t, float u, float w, float v, float s) {
    EigenType::Matrix6f c;
    c << u, w, -w, x, s, s,
        -t, y, v, -v, z, -t,
        t, v, y, v, t, -z,
        z, -v, v, y, -t, -t,
        s, x, -w, w, u, s,
        s, w, -x, w, s, u;
    return c;
}

bool DihedralGroupOnSix::checkIfInCentralizer(EigenType::Matrix6f m) {
    return m.isApprox(matrixFormOfCentralizer(m(0,3), m(1,1), m(3, 0), m(2, 0), m(0, 0), m(0, 1), m(1, 2), m(0, 4)));
}

std::string DihedralGroupOnSix::groupName() {
    return "DihedralGroupD_6";
}
