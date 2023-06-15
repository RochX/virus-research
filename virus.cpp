//
// Created by xavier on 6/13/2023.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <omp.h>

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#define NUM_THREADS 16

using namespace std;

typedef Eigen::Vector<float, 6> Vector6f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;

void generateAllCentralizerCandidates(const string&, const string&, bool);
bool readnextMatrix(ifstream&, Matrix6f&);
bool readnextNMatrices(ifstream&, vector<Matrix6f>&, int);
vector<Matrix6f> convertLinesToMatrices(const vector<string>&);
bool readnextNLines(ifstream&, vector<string>&, int);
void fileOutputAllFullRankMatrices(const std::vector<vector<Vector6f>>&, const string&, bool);
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
            break;
        }
    }

    Matrix6f D;
    D << 1,1,-1,-1,1,1,
        1,1,1,-1,-1,1,
        -1,1,1,1,-1,1,
        -1,-1,1,1,1,1,
        1,-1,-1,1,1,1,
        1,1,1,1,1,1;
    D *= 0.5;

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

    std::vector<Vector6f> orbit_s, orbit_b, orbit_Dinvs, orbit_Dinvb, orbit_Dinvf, P_0, P_1;
    orbit_s = generateOrbit(s, ICO);
    orbit_b = generateOrbit(b, ICO);
    orbit_Dinvb = generateOrbit(D.inverse()*b, ICO);
    orbit_Dinvf = generateOrbit(D.inverse()*f, ICO);

    // create P_0 by appending all the orbits together
    append_vector(P_0, orbit_s);
    append_vector(P_0, orbit_b);
    append_vector(P_0, orbit_Dinvb);
    append_vector(P_0, orbit_Dinvf);

    // create the list that holds all the possibilities for each column vector
    std::vector<std::vector<Vector6f>> B0_choices;
    B0_choices.push_back(orbit_s);
    B0_choices.push_back(orbit_Dinvb);
    B0_choices.push_back(orbit_Dinvf);
    B0_choices.push_back(orbit_b);
    B0_choices.push_back(P_0);
    B0_choices.push_back(P_0);


    fileOutputAllFullRankMatrices(B0_choices, "B0_matrices.csv", true);

    // process of picking linearly independent matrices from P_1
    orbit_s = generateOrbit(s, ICO);
    orbit_b = generateOrbit(b, ICO);
    orbit_Dinvs = generateOrbit(D.inverse()*s, ICO);
    orbit_Dinvb = generateOrbit(D.inverse()*b, ICO);
    orbit_Dinvf = generateOrbit(D.inverse()*f, ICO);

    append_vector(P_1, orbit_b);
    append_vector(P_1, orbit_Dinvs);
    append_vector(P_1, orbit_Dinvb);
    append_vector(P_1, orbit_Dinvf);
    append_vector(P_1, orbit_s);

    // append choices for B1
    std::vector<vector<Vector6f>> B1_choices;
    B1_choices.push_back(orbit_b);
    B1_choices.push_back(orbit_Dinvs);
    B1_choices.push_back(orbit_Dinvb);
    B1_choices.push_back(orbit_Dinvf);
    B1_choices.push_back(orbit_s);
    B1_choices.push_back(P_1);

    fileOutputAllFullRankMatrices(B1_choices, "B1_matrices.csv", true);
}

// using two filenames, generate B_1B_0^-1 and check if it's in the centralizer
void generateAllCentralizerCandidates(const string &B0_filename, const string &B1_filename, bool permuteB1) {
    ifstream b0in (B0_filename);
    ifstream b1in (B1_filename);

    string b0line, b1line;
    string cell;
    stringstream lineStream;
    vector<float> entries;

    // read csv headers
    getline(b0in, b0line);
    getline(b1in, b1line);

    vector<Matrix6f> b1matrices;
    Matrix6f b0matrix;

    int tempCount = 0;

    int b0count = 0;
    int b1count = 0;
    int b0CAP = 10;
    int b1CAP = 100000;
    int N = 10000;
    int ICO_centralizer_count = 0;

    // read B0 matrices
    while (readnextMatrix(b0in, b0matrix)) {

        // reset stuff for B1
        b1count = 0;
        b1in.clear();
        b1in.seekg(0);
        // read csv header
        getline(b1in, b1line);

        // read B1 matrices
        while(readnextNMatrices(b1in, b1matrices, N)) {
            // check for ICO centralizer using B_1B_0^-1
            for (const Matrix6f &b1 : b1matrices) {
//                cout << b1*b0matrix.inverse() << endl << endl;
                if (check_ICO_centralizer(b1*b0matrix.inverse())) {
                    ICO_centralizer_count++;
                }
            }


            // keep track of how many B1 matrices read
            b1count += N;
            if (b1count >= b1CAP)
                break;
        }

        // keep track of how many B0 matrices read
        b0count++;
        cout << b0count << endl;
        if (b0count >= b0CAP)
            break;
    }

    cout << "ICO centralizer count:\t" << ICO_centralizer_count << endl;

    b0in.close();
    b1in.close();
}

bool readnextMatrix(ifstream &fin, Matrix6f &m) {
    vector<Matrix6f> matrices;
    // try to read matrix
    if (readnextNMatrices(fin, matrices, 1)) {
        // succeeded, assign it
        m = matrices.front();
        return true;
    }
    // failed to read matrix
    else {
        return false;
    }
}

bool readnextNMatrices(ifstream &fin, vector<Matrix6f> &matrices, int N) {
    vector<string> lines;

    readnextNLines(fin, lines, N);

    matrices = convertLinesToMatrices(lines);

    return !matrices.empty();
}

// takes a list of lines and converts them to a list of 6x6 matrices
// precondition: each string in the vector lines is 36 comma delimited floats
vector<Matrix6f> convertLinesToMatrices(const vector<string> &lines) {
    vector<Matrix6f> matrices;
    Matrix6f m;

    stringstream ss;
    string cell;
    int pos;

    for (const string &line : lines) {
        ss.clear();
        ss << line;

        pos = 0;
        while(getline(ss, cell, ',')) {
            m(pos/6, pos%6) = stof(cell);
            pos++;
        }

        matrices.push_back(m);
    }

    return matrices;
}

// reads the next N lines from the specified ifstream
// the ifstream has less than N lines left, it will read the rest of them
bool readnextNLines(ifstream &fin, vector<string> &lines, int N) {
    string line;
    int lines_read = 0;
    while (getline(fin, line)) {
        lines.push_back(line);

        lines_read++;
        if (lines_read >= N) {
            break;
        }
    }

    return lines_read > 0;
}

// try all possibilities of making 6x6 matrices using our set P
// precondition: P has precisely 6 lists of 6 element column vectors
// precondition: filename is a csv file
void fileOutputAllFullRankMatrices(const std::vector<vector<Vector6f>> &P, const string &filename, bool parallelize) {
    string CSV_HEADER = "11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 61, 62, 63, 64, 65, 66";
    Eigen::IOFormat CommaSepVals(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");

    Matrix6f m;
    int threads = (parallelize) ? NUM_THREADS : 1;
    long long currChecked;
    long long totalPmatrices = 1;
    long long totalrank6Pmatrices = 0;

    for (const vector<Vector6f> &v : P) {
        totalPmatrices *= v.size();
    }

    // check we actually have a nonzero number of matrices to check
    if (totalPmatrices < 1) {
        cerr << "0 P matrices! One of the lists within P has length zero!" << endl;
        return;
    }

    // check if we opened file properly, if so output csv header
    ofstream fout (filename);
    if (!fout.is_open()) {
        cerr << "Failed to open file " << filename << " within function fileOutputAllFullRankMatrices." << endl;
        return;
    }
    else {
        fout << CSV_HEADER << endl;
    }

    // start tracking time
    auto start_time = omp_get_wtime();

    // do parallel computation
    #pragma omp parallel num_threads(threads) private(m, currChecked) shared(CommaSepVals, fout, P, totalPmatrices, totalrank6Pmatrices, std::cout) default(none)
    {
        currChecked = 0;

        // brute force approach
        #pragma omp for collapse(6) reduction(+:totalrank6Pmatrices)
        for (int first = 0; first < P[0].size(); ++first) {
            for (int second = 0; second < P[1].size(); ++second) {
                for (int third = 0; third < P[2].size(); ++third) {
                    for (int fourth = 0; fourth < P[3].size(); ++fourth) {
                        for (int fifth = 0; fifth < P[4].size(); ++fifth) {
                            for (int sixth = 0; sixth < P[5].size(); ++sixth) {
                                m << P[0][first], P[1][second], P[2][third], P[3][fourth], P[4][fifth], P[5][sixth];

                                // check our matrix is of rank 6
                                if (m.colPivHouseholderQr().rank() == 6) {
                                    #pragma omp critical
                                    fout << m.format(CommaSepVals) << endl;

                                    // count how many matrices we output
                                    totalrank6Pmatrices++;
                                }

                                // get a sense of progress done using thread 0
                                if (omp_get_thread_num() == 0) {
                                    currChecked++;
                                    if (currChecked % 1000 == 0) {
                                        cout << "Thread 0 has checked:\t" << currChecked << " out of " << totalPmatrices/omp_get_num_threads() << endl;
                                    }
                                    else if (currChecked == totalPmatrices/omp_get_num_threads()) {
                                        cout << "Thread 0 has checked:\t" << currChecked << " out of " << totalPmatrices/omp_get_num_threads() << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // finished with parallel region means we finished with file
    fout.close();

    // finish tracking time, output some results
    auto current_time = omp_get_wtime();
    cout << "\nFinished outputting to " << filename << "." << endl;
    cout << "Outputted " << totalrank6Pmatrices << " matrices of rank 6." << endl;
    cout << "Total time taken:\t" << (current_time - start_time) << " seconds" << endl;
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