#include "DihedralGroupOnTen.hpp"

DihedralGroupOnTen::DihedralGroupOnTen() {
    // A has order 5
    generatorA << 0, 0, 0, 0, 1, 0,
            1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;

    // B has order 2
    generatorB << 0, -1, 0, 0, 0, 0,
            -1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, -1, 0,
            0, 0, 0, -1, 0, 0,
            0, 0, -1, 0, 0, 0,
            0, 0, 0, 0, 0, -1;

    // D_10 has 10 elements
    createGroupFromTwoGeneratorsAndKnownSize(generatorA, generatorB, 10);
}

EigenType::Matrix6f DihedralGroupOnTen::matrixFormOfCentralizer(float x, float y, float z, float t, float u, float w) {
    EigenType::Matrix6f c;
    c << z, x, y, y, x, t,
        x, z, x, y, y, t,
        y, x, z, x, y, t,
        y, y, x, z, x, t,
        x, y, y, x, z, t,
        u, u, u, u, u, w;
    return c;
}

bool DihedralGroupOnTen::checkIfInCentralizer(EigenType::Matrix6f m) {
    return m.isApprox(matrixFormOfCentralizer(m(0,1), m(0,2), m(0,0), m(0, 5), m(5,0), m(5,5)));
}

std::string DihedralGroupOnTen::groupName() {
    return "DihedralGroupD_10";
}
