//
// Created by xavier on 6/13/2023.
//
#include <cstdio>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <set>

using namespace std;

typedef Eigen::Vector<float, 6> Vector6f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;

std::vector<Vector6f> generateOrbit(const Vector6f&, std::vector<Matrix6f>&);

int main() {
    // initialize the generator matrices of the group ICO
    Matrix6f A, B, ID;
    ID = Matrix6f::Identity();

    // A has order 2
    A << -1, 0, 0, 0, 0, 0,
        0, -1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0;

    // B has order 3
    B << 0, -1, 0, 0, 0, 0,
        0, 0, -1, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0;

    std::cout << "A:\n" << A << std::endl << std::endl;
    std::cout << "B:\n" << B << std::endl << std::endl;

    // check order of matrices;
    std::cout << "A^2:\n" << A*A << std::endl << std::endl;
    std::cout << "B^3:\n" << B*B*B << std::endl << std::endl;

    // note that this "vector" is a C++ dynamic list, *not* the Eigen library vector.
    std::vector<Matrix6f> ICO;
    Matrix6f curr;
    ICO.push_back(ID);
    int element_id;
    // generate ICO with a brute force approach as we know it has 60 elements.
    for (int i = 0; i < 10000; ++i) {
        element_id = i;
        curr = ID;

        // generate ICO by using binary representation, with 0 corresponding to A and 1 corresponding to B.
        // exs: 10 = BA
        //      100 = BA^2
        //      1011 = BAB^2
        // note that while we do not get elements like 011 = AB^2 explicitly, we know that B has order 3, and so we do get the element 111011 = B^3AB^2 = AB^2
        // so this process will generate all elements of ICO.
        do {
            if (element_id % 2 == 0) {
                curr *= A;
            }
            else {
                curr *= B;
            }
            element_id /= 2;
        } while (element_id > 0);

        // check if ICO already has the current element, if not add it.
        if (std::find(ICO.begin(), ICO.end(), curr) == ICO.end())
            ICO.push_back(curr);

        // if ICO has 60 (unique) elements we're done
        if (ICO.size() == 60) {
            cout << i << endl;
            break;
        }
    }

//    cout << "ICO:\n";
//    for (Matrix6f m : ICO) {
//        //cout << m << endl << endl;
//        // outputs the matrix in a 1D view, row wise.
//        cout << m.reshaped<Eigen::RowMajor>().transpose() << endl;
//    }
//    cout << ICO.size() << endl;

    Matrix6f D;
    D << 1,1,-1,-1,1,1,
        1,1,1,-1,-1,1,
        -1,1,1,1,-1,1,
        -1,-1,1,1,1,1,
        1,-1,-1,1,1,1,
        1,1,1,1,1,1;
    D *= 0.5;
//    cout << D << endl;

    // declare and initialize vectors which we will make orbits of
    // s: the ICO orbit of s is the icosahedron (the particular 6D embedding we are working with)
    // b: the ICO orbit of b is the dodecahedron (which is dual to the icosahedron which is why this vector by itself produces an SC lattice -- in fact, the specific SC lattice is D*sc where sc denotes our standard(fixed) simple cubic)
    // f: the ICO orbit of f is the icosidodecahedron
    Vector6f s, b, f;
    s << 1, 0, 0, 0, 0, 0;
    b << 1,-1,1,1,-1,1;
    b *= 0.5;
    f << 1,0,0,-1,0,0;
    f *= 0.5;

    std::vector<Vector6f> orbit_s, orbit_b, orbit_Dinvb, orbit_Dinvf;
    orbit_s = generateOrbit(s, ICO);
    orbit_b = generateOrbit(b, ICO);
    orbit_Dinvb = generateOrbit(D.inverse()*b, ICO);
    orbit_Dinvf = generateOrbit(D.inverse()*f, ICO);

    for (Vector6f v : orbit_s) {
        cout << v.transpose() << endl;
    }
}

// function to generate orbit of a vector under a given group.
std::vector<Vector6f> generateOrbit(const Vector6f &v, std::vector<Matrix6f> &G) {
    std::vector<Vector6f> orbit;

    for (const Matrix6f &g : G) {
        orbit.emplace_back(g*v);
    }

    return orbit;
}