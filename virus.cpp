//
// Created by xavier on 6/13/2023.
//
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <chrono>

using namespace std;

typedef Eigen::Vector<float, 6> Vector6f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;

void append_vector(std::vector<Vector6f>&, std::vector<Vector6f>&, bool = true);
std::vector<Vector6f> generateOrbit(const Vector6f&, std::vector<Matrix6f>&);
Matrix6f ICO_centralizer(float, float);
bool check_ICO_centralizer(Matrix6f);

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

    std::vector<Vector6f> orbit_s, orbit_b, orbit_Dinvb, orbit_Dinvf, P_0;
    orbit_s = generateOrbit(s, ICO);
    orbit_b = generateOrbit(b, ICO);
    orbit_Dinvb = generateOrbit(D.inverse()*b, ICO);
    orbit_Dinvf = generateOrbit(D.inverse()*f, ICO);

    // create P_0 by appending all the orbits together
    append_vector(P_0, orbit_s);
    append_vector(P_0, orbit_b);
    append_vector(P_0, orbit_Dinvb);
    append_vector(P_0, orbit_Dinvf);

    for (Vector6f v : orbit_s) {
        cout << v.transpose() << endl;
    }

    // TODO: add the picking of linearly independent matrices from orbits and P_0
    std::vector<std::vector<Vector6f>> B0_choices;
    B0_choices.push_back(orbit_s);
    B0_choices.push_back(orbit_Dinvb);
    B0_choices.push_back(orbit_Dinvf);
    B0_choices.push_back(orbit_b);
    B0_choices.push_back(P_0);
    B0_choices.push_back(P_0);





    long long totalB0matrices = 1;
    for (const vector<Vector6f>& mv : B0_choices) {
        totalB0matrices *= mv.size();
    }
    cout << totalB0matrices << endl;

    std::vector<Matrix6f> rank6matrices;
    Matrix6f m;
    m.setZero();
    int rank;
    // this part sucks because 6 nested for loops, but it's the simplest way to approach this
    auto start_time = std::chrono::high_resolution_clock::now();
    cout << "Starting nested loops...\n";
    for (int first = 0; first < B0_choices[0].size(); ++first) {
        for (int second = 0; second < B0_choices[1].size(); ++second) {
            for (int third = 0; third < B0_choices[2].size(); ++third) {
                for (int fourth = 0; fourth < B0_choices[3].size(); ++fourth) {
                    for (int fifth = 0; fifth < B0_choices[4].size(); ++fifth) {
                        for (int sixth = 0; sixth < B0_choices[5].size(); ++sixth) {
                            m << B0_choices[0][first], B0_choices[1][second], B0_choices[2][third], B0_choices[3][fourth], B0_choices[4][fifth], B0_choices[5][sixth];

                            // check our matrix is of rank 6
                            if (m.colPivHouseholderQr().rank() == 6) {
                                rank6matrices.push_back(m);
//                                if (rank6matrices.size() >= 100) {
//                                    cout << first << "\t" << second << "\t" << third << "\t" << fourth << "\t" << fifth << "\t" << sixth << endl;
//                                    goto doneloop;
//                                }
                            }
                        }
                    }
                }
            }
        }
    }
//    doneloop:
    cout << "\nDone.\n";
    auto current_time = std::chrono::high_resolution_clock::now();

    std::cout << "Loop ran for: " << (current_time - start_time).count()/1000000000.0 << " seconds" << std::endl;

//    for (Matrix6f M : rank6matrices) {
//        cout << M << endl << endl;
//    }
    cout << rank6matrices.size() << endl;




    // TODO: add the picking of linearly independent matrices from orbits and P_1 (and make P_1)
}

// append the contents of v2 to v1, possibly checking for duplicates
void append_vector(std::vector<Vector6f>& v1, std::vector<Vector6f>& v2, bool removeDups) {
    for (const Vector6f &v : v2) {
        // if bool is true only append if it is not a duplicate
        if (removeDups) {
            if (std::find(v1.begin(), v1.end(), v) == v1.end()) {
                v1.push_back(v);
            }
        }
        // otherwise just append without checking
        else {
            v1.push_back(v);
        }
    }
}

// function to generate orbit of a vector under a given group.
std::vector<Vector6f> generateOrbit(const Vector6f &v, std::vector<Matrix6f> &G) {
    std::vector<Vector6f> orbit;

    for (const Matrix6f &g : G) {
        // ensure no duplicates in orbit
        if (std::find(orbit.begin(), orbit.end(), g*v) == orbit.end())
            orbit.emplace_back(g*v);
    }

    return orbit;
}

Matrix6f ICO_centralizer(float z, float x) {
    Matrix6f c;
    c << z,x,-x,-x,x,x,
        x,z,x,-x,-x,x,
        -x,x,z,x,-x,x,
        -x,-x,x,z,x,x,
        x,-x,-x,x,z,x,
        x,x,x,x,x,z;
    return c;
}

bool check_ICO_centralizer(Matrix6f m) {
    return m == ICO_centralizer(m(0,0), m(0,1));
}