#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <omp.h>

#include "EigenTypes.hpp"
#include "IcosahedralGroup.hpp"
#include "Matrix6fFileReader.hpp"
#include "TetrahedralGroup.hpp"

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#define NUM_THREADS 16

using namespace EigenType;

void generateAllCentralizerCandidates(Matrix6fGroup&, const std::string&, const std::string&, bool);
void generateAllB0andB1Matrices(Matrix6fGroup&, bool, bool, bool = false);
void fileOutputAllFullRankMatrices(const std::vector<std::vector<Vector6f>>&, const std::string&, bool);
void append_vector(std::vector<Vector6f>&, std::vector<Vector6f>&, bool = true);

std::string B0MatricesFileName;
std::string B1MatricesFileName;

int main(int argc, char *argv[]) {
    std::string currDirectory;
    // if second argument given, assume all files are in directory specified by the second argument
    if (argc > 1) {
        currDirectory = argv[1];
        if (currDirectory.back() != '/') {
            currDirectory = currDirectory + "/";
        }
    }

    TetrahedralGroup TetrahedralGroup;
    IcosahedralGroup IcosahedralGroup;

    Matrix6fGroup* currentGroup = &TetrahedralGroup;

    B0MatricesFileName = currDirectory + currentGroup->groupName() + "_B0_matrices.csv";
    B1MatricesFileName = currDirectory + currentGroup->groupName() + "_B1_matrices.csv";

    generateAllB0andB1Matrices(*currentGroup, true, true, true);

    generateAllCentralizerCandidates(*currentGroup, B0MatricesFileName, B1MatricesFileName, true);
}

// using two filenames, generate B_1B_0^-1 and check if it's in the centralizer
void generateAllCentralizerCandidates(Matrix6fGroup &matrixGroup, const std::string &B0_filename, const std::string &B1_filename, bool permuteB1) {
    std::ifstream b0in (B0_filename);
    std::ifstream b1in (B1_filename);

    std::string b0line, b1line;
    std::string cell;
    std::stringstream lineStream;
    std::vector<float> entries;

    // read csv headers
    getline(b0in, b0line);
    getline(b1in, b1line);

    std::vector<Matrix6f> b1matrices;
    Matrix6f b0matrix;

    int tempCount = 0;

    int b0count = 0;
    int b1count = 0;
    int b0CAP = 10;
    int b1CAP = 100000;
    int N = 10000;
    int group_centralizer_count = 0;

    // read B0 matrices
    while (Matrix6fFileReader::readNextMatrix(b0in, b0matrix)) {

        // reset stuff for B1
        b1count = 0;
        b1in.clear();
        b1in.seekg(0);
        // read csv header
        getline(b1in, b1line);

        // read B1 matrices
        while(Matrix6fFileReader::readNextNMatrices(b1in, b1matrices, N)) {
            // check for ICO centralizer using B_1B_0^-1
            for (const Matrix6f &b1 : b1matrices) {
//                std::cout << b1*b0matrix.inverse() << std::endl << std::endl;
                if (matrixGroup.checkIfInCentralizer(b1*b0matrix.inverse())) {
                    group_centralizer_count++;
                }
            }


            // keep track of how many B1 matrices read
            b1count += N;
            if (b1count >= b1CAP)
                break;
        }

        // keep track of how many B0 matrices read
        b0count++;
        std::cout << b0count << std::endl;
        if (b0count >= b0CAP)
            break;
    }

    std::cout << matrixGroup.groupName() << " centralizer count:\t" << group_centralizer_count << std::endl;

    b0in.close();
    b1in.close();
}

void generateAllB0andB1Matrices(Matrix6fGroup &matrixGroup, bool generateB0, bool generateB1, bool generatePartial) {
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
    orbit_s = matrixGroup.getOrbitOfVector(s);
    orbit_b = matrixGroup.getOrbitOfVector(b);
    orbit_Dinvb = matrixGroup.getOrbitOfVector(D.inverse() * b);
    orbit_Dinvf = matrixGroup.getOrbitOfVector(D.inverse() * f);

    // create P_0 by appending all the orbits together
    append_vector(P_0, orbit_s);
    append_vector(P_0, orbit_b);
    append_vector(P_0, orbit_Dinvb);
    append_vector(P_0, orbit_Dinvf);

    if (generatePartial) {
        P_0.clear();
        P_0.push_back(s);
        P_0.push_back(b);
    }

    // create the list that holds all the possibilities for each column vector
    std::vector<std::vector<Vector6f>> B0_choices;
    B0_choices.push_back(orbit_s);
    B0_choices.push_back(orbit_Dinvb);
    B0_choices.push_back(orbit_Dinvf);
    B0_choices.push_back(orbit_b);
    B0_choices.push_back(P_0);
    B0_choices.push_back(P_0);



    // process of picking linearly independent matrices from P_1
    orbit_s = matrixGroup.getOrbitOfVector(s);
    orbit_b = matrixGroup.getOrbitOfVector(b);
    orbit_Dinvs = matrixGroup.getOrbitOfVector(D.inverse()*s);
    orbit_Dinvb = matrixGroup.getOrbitOfVector(D.inverse()*b);
    orbit_Dinvf = matrixGroup.getOrbitOfVector(D.inverse()*f);

    append_vector(P_1, orbit_b);
    append_vector(P_1, orbit_Dinvs);
    append_vector(P_1, orbit_Dinvb);
    append_vector(P_1, orbit_Dinvf);
    append_vector(P_1, orbit_s);

    if (generatePartial) {
        P_1.clear();
        P_1.push_back(s);
        P_1.push_back(b);
    }

    // append choices for B1
    std::vector<std::vector<Vector6f>> B1_choices;
    B1_choices.push_back(orbit_b);
    B1_choices.push_back(orbit_Dinvs);
    B1_choices.push_back(orbit_Dinvb);
    B1_choices.push_back(orbit_Dinvf);
    B1_choices.push_back(orbit_s);
    B1_choices.push_back(P_1);

    if (generateB0)
        fileOutputAllFullRankMatrices(B0_choices, B0MatricesFileName, true);

    if (generateB1)
        fileOutputAllFullRankMatrices(B1_choices, B1MatricesFileName, true);
}


// try all possibilities of making 6x6 matrices using our set P
// precondition: P has precisely 6 lists of 6 element column vectors
// precondition: filename is a csv file
void fileOutputAllFullRankMatrices(const std::vector<std::vector<Vector6f>> &P, const std::string &filename, bool parallelize) {
    std::string CSV_HEADER = "11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 61, 62, 63, 64, 65, 66";
    Eigen::IOFormat CommaSepVals(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");

    Matrix6f m;
    int threads = (parallelize) ? NUM_THREADS : 1;
    long long currChecked;
    long long totalPmatrices = 1;
    long long totalrank6Pmatrices = 0;

    for (const std::vector<Vector6f> &v : P) {
        totalPmatrices *= v.size();
    }

    // check we actually have a nonzero number of matrices to check
    if (totalPmatrices < 1) {
        std::cerr << "0 P matrices! One of the lists within P has length zero!" << std::endl;
        return;
    }

    // check if we opened file properly, if so output csv header
    std::ofstream fout (filename);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file " << filename << " within function fileOutputAllFullRankMatrices." << std::endl;
        return;
    }
    else {
        fout << CSV_HEADER << std::endl;
    }

    std::cout << "Starting full rank matrix output to file:\t" + filename << std::endl;

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
                                    fout << m.format(CommaSepVals) << std::endl;

                                    // count how many matrices we output
                                    totalrank6Pmatrices++;
                                }

                                // get a sense of progress done using thread 0
                                if (omp_get_thread_num() == 0) {
                                    currChecked++;
                                    if (currChecked % 1000 == 0) {
                                        std::cout << "Thread 0 has checked:\t" << currChecked << " out of " << totalPmatrices/omp_get_num_threads() << std::endl;
                                    }
                                    else if (currChecked == totalPmatrices/omp_get_num_threads()) {
                                        std::cout << "Thread 0 has checked:\t" << currChecked << " out of " << totalPmatrices/omp_get_num_threads() << std::endl;
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
    std::cout << "\nFinished outputting to " << filename << "." << std::endl;
    std::cout << "Outputted " << totalrank6Pmatrices << " matrices of rank 6." << std::endl;
    std::cout << "Total time taken:\t" << (current_time - start_time) << " seconds" << std::endl << std::endl;
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