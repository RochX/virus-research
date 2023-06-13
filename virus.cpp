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

typedef Eigen::Matrix<float, 6, 6> Matrix6f;

int main() {
    // initialize the generator matrices of the group ICO
    Matrix6f A, B, ID;
    ID << 1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1;

    A << -1, 0, 0, 0, 0, 0,
        0, -1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0;

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
    vector<Matrix6f> ICO;
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

    cout << "ICO:\n";
    for (Matrix6f m : ICO) {
        cout << m << endl << endl;
    }
    cout << ICO.size() << endl;

}