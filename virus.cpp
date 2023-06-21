#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <omp.h>

#include "EigenTypes.hpp"
#include "IcosahedralGroup.hpp"
#include "Matrix6fFileReader.hpp"
#include "PermutationGroup.hpp"
#include "TetrahedralGroup.hpp"
#include "std_vector_functions.hpp"

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#define NUM_THREADS 16

using namespace EigenType;

std::vector<Matrix6f> GetAllPossibleTMatricesInIco(const std::vector<float>& T_entries);
void generateAllCentralizerCandidates(Matrix6fGroup&, const std::string&, const std::string&, bool);
void generateAllB0andB1Matrices(Matrix6fGroup&, bool, bool, bool = false);
void fileOutputAllFullRankMatrices(const std::vector<std::vector<Vector6f>>&, const std::string&, bool);

// TODO: move these functions into their own file...
std::vector<float> findPossibleTEntriesWithOneSample(const std::string &B0_filename, const std::string &B1_filename, int sample_size, int skip);
std::vector<float> findPossibleTEntriesBySampling(const std::string &B0_filename, const std::string &B1_filename, int num_samples, int sample_size);
std::vector<float> findPossibleTEntriesHardCoded();
std::vector<Matrix6f> GetAllPossibleTMatricesInIco(const std::vector<float>& T_entries);

bool AllEntriesOfParticularValues(const Matrix6f& matrix, const std::vector<float>& valid_values);
bool FloatsAreApproxEqual(float x, float y);

std::string B0MatricesFileName;
std::string B1MatricesFileName;
std::string CentralizerFileName;

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

    // always generate B0 and B1 matrices using ICO
    Matrix6fGroup* currentGroup = &IcosahedralGroup;
    B0MatricesFileName = currDirectory + currentGroup->groupName() + "_B0_matrices.csv";
    B1MatricesFileName = currDirectory + currentGroup->groupName() + "_B1_matrices.csv";
    generateAllB0andB1Matrices(*currentGroup, true, true, true);

//    // here the group may vary
//    currentGroup = &TetrahedralGroup;
//    CentralizerFileName = currDirectory + currentGroup->groupName() + "_centralizer_matrices.csv";
//    generateAllCentralizerCandidates(*currentGroup, B0MatricesFileName, B1MatricesFileName, true);

    std::vector<float> T_entries;
//    T_entries = findPossibleTEntriesBySampling(B0MatricesFileName, B1MatricesFileName, 10, 1000);
    T_entries = findPossibleTEntriesHardCoded();

    std::vector<Matrix6f> T_matrices = GetAllPossibleTMatricesInIco(T_entries);

    std::ifstream b0in (B0MatricesFileName);
    std::string csv_header;
    getline(b0in, csv_header);

    std::vector<float> B1_vals {-1, -0.5, 0, 0.5, 1};
    long long invalid_count = 0;
    Matrix6f B0_m;
    const int NUM_CHECK = 10;
    for (int i = 0; i < NUM_CHECK && Matrix6fFileReader::readNextMatrix(b0in, B0_m); i++) {
        for (const Matrix6f& T : T_matrices) {
            if (!AllEntriesOfParticularValues(T*B0_m, B1_vals)) {
                invalid_count++;
            }
        }
    }
    std::cout << "Number invalid is: " << invalid_count << " out of " << T_matrices.size()*NUM_CHECK << std::endl;
}

std::vector<Matrix6f> GetAllPossibleTMatricesInIco(const std::vector<float>& T_entries) {
    std::vector<Matrix6f> T_matrices_in_ICO;
    for (float f1 : T_entries) {
        for (float f2: T_entries) {
            T_matrices_in_ICO.push_back(IcosahedralGroup::matrixFormOfCentralizer(f1, f2));
        }
    }

    return T_matrices_in_ICO;
}

std::vector<float> findPossibleTEntriesHardCoded() {
    std::vector<float> possible_T_entries;

    float sixths = -10;
    while (sixths <= 10) {
        std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, sixths, 0.0001);
        sixths += 1.0f/6;
    }

    float fourths = -10;
    while (fourths <= 10) {
        std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, fourths, 0.0001);
        fourths += 1.0f/4;
    }

    std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, 10, 0.0001);

    std::sort(possible_T_entries.begin(), possible_T_entries.end());
    return possible_T_entries;
}

std::vector<float> findPossibleTEntriesBySampling(const std::string &B0_filename, const std::string &B1_filename, int num_samples, int sample_size) {
    std::vector<float> possible_T_entries;

    for (int k = 0; k < num_samples; ++k) {
        std::cout << "Sample " << k << ": ";
        std::vector<float> T_entries = findPossibleTEntriesWithOneSample(B0MatricesFileName, B1MatricesFileName,
                                                                         sample_size, k * sample_size);
        for (float f: T_entries) {
            std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, f, 0.0001);
        }
    }

    std::sort(possible_T_entries.begin(), possible_T_entries.end());
    return possible_T_entries;
}

std::vector<float> findPossibleTEntriesWithOneSample(const std::string &B0_filename, const std::string &B1_filename, int sample_size, int skip) {
    std::cout << "Finding all possible T entries..." << std::endl;
    auto start_time = omp_get_wtime();

    std::ifstream b0in (B0_filename);
    std::ifstream b1in (B1_filename);

    std::string csv_header;
    getline(b0in, csv_header);
    getline(b1in, csv_header);

    Matrix6f dump;
    // for simplicity just check a sample_size * sample_size, but with an option to skip through the file.
    for (int i = 0; i < skip; ++i) {
        Matrix6fFileReader::readNextMatrix(b0in, dump);
        Matrix6fFileReader::readNextMatrix(b1in, dump);
    }

    // check if sample_size is too big, and cap it if yes with warning...
    int MAX_SAMPLE_SIZE = 1000000;
    if (sample_size >= MAX_SAMPLE_SIZE) {
        std::cout << "WARNING: Sample size over 1 million, capping to 1 million." << std::endl;
        sample_size = MAX_SAMPLE_SIZE;
    }

    std::vector<float> possible_T_entries;
    std::vector<Matrix6f> B0_matrices_list, B1_matrices_list;
    Matrix6fFileReader::readNextNMatrices(b0in, B0_matrices_list, sample_size);
    Matrix6fFileReader::readNextNMatrices(b1in, B1_matrices_list, sample_size);
    #pragma omp parallel num_threads(NUM_THREADS) shared(possible_T_entries, B0_matrices_list, B1_matrices_list) default(none)
    {
        const float EPSILON = 0.0001;
        std::vector<float> thread_possible_T_entries;
        #pragma omp for collapse(2)
        for (int i = 0; i < B0_matrices_list.size(); ++i) {
            for (int j = 0; j < B1_matrices_list.size(); ++j) {
                Matrix6f B0_m = B0_matrices_list[i];
                Matrix6f B1_m = B1_matrices_list[j];
                Matrix6f T = B1_m * B0_m.inverse();
                bool add_current_entry;

                // check if any of the entries of T are already accounted for.
                for (int k = 0; k < 36; ++k) {
                    std_vector_functions::push_backIfNotInVector<float>(thread_possible_T_entries, T(k / 6, k % 6), 0.0001);
                }
            }
        }

        // combine into one list
        #pragma omp critical
        {
            bool add_to_list;
            for (float f : thread_possible_T_entries) {
                add_to_list = true;
                for (float t : possible_T_entries) {
                    if (std::fabs(f - t) < EPSILON) {
                        add_to_list = false;
                        break;
                    }
                }

                if (add_to_list)
                    possible_T_entries.push_back(f);
            }
        }
    }

    std::sort(possible_T_entries.begin(), possible_T_entries.end());
    auto current_time = omp_get_wtime();
    std::cout << "Done finding T entries, done in " << (current_time-start_time) << std::endl;

    return possible_T_entries;
}

// using two filenames, generate B_1B_0^-1 and check if it's in the centralizer
void generateAllCentralizerCandidates(Matrix6fGroup &matrixGroup, const std::string &B0_filename, const std::string &B1_filename, bool permuteB1) {
    PermutationGroup permutationGroup;

    std::ifstream b0in (B0_filename);
    std::ifstream b1in (B1_filename);
    std::ofstream centralizerFout(CentralizerFileName);
    Eigen::IOFormat CommaSepVals(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");

    std::string b0line, b1line;
    std::string cell;
    std::stringstream lineStream;

    // read csv headers
    getline(b0in, b0line);
    getline(b1in, b1line);

    std::vector<Matrix6f> b1matrices;
    Matrix6f b0matrix;
    std::vector<Matrix6f> permutationMatrices = permutationGroup.getGroupElements();

    if (!permuteB1) {
        permutationMatrices.clear();
        permutationMatrices.emplace_back(Matrix6f::Identity());
    }

    int b0count, b1count;
    int b0CAP = INT32_MAX;
    int b1CAP = INT32_MAX;
    int N = 10000;
    int group_centralizer_count = 0;

    std::cout << "Start checking all pairs of B0 and B1..." << std::endl;
    auto start_time = omp_get_wtime();

    #pragma omp parallel num_threads(NUM_THREADS) private(b0count, b1count, b0matrix, b1matrices, b1line) shared(permutationMatrices, CommaSepVals, centralizerFout, matrixGroup, start_time, std::cout, N, b0CAP, b1CAP, b0in, b1in, group_centralizer_count) default(none)
    {
        #pragma omp master
        {
            b0count = 0;
            // read B0 matrices
            while (Matrix6fFileReader::readNextMatrix(b0in, b0matrix)) {

                // reset B1 file reader to top of file
                b1count = 0;
                b1in.clear();
                b1in.seekg(0);
                // read csv header
                getline(b1in, b1line);

                // read B1 matrices
                while (Matrix6fFileReader::readNextNMatrices(b1in, b1matrices, N)) {
                    // check for ICO centralizer using B_1B_0^-1
                    for (const Matrix6f &b1: b1matrices) {
                        // assign task to check this pair of matrices (with permutation)
                        #pragma omp task shared(std::cout, permutationMatrices, centralizerFout, CommaSepVals, matrixGroup, group_centralizer_count) firstprivate(b0matrix, b1) default(none)
                        {
                            for (const Matrix6f &P: permutationMatrices) {
                                {
                                    if (matrixGroup.checkIfInCentralizer(b1 * P * b0matrix.inverse())) {
                                        #pragma omp critical
                                        {
                                            centralizerFout << b1.format(CommaSepVals) << std::endl
                                                            << b0matrix.format(CommaSepVals) << std::endl << std::endl;
                                            group_centralizer_count++;
                                        }
                                    }
                                }
                            }
                        }
                    }


                    // keep track of how many B1 matrices read
                    b1count += N;
                    if (b1count >= b1CAP)
                        break;
                }

                // keep track of how many B0 matrices read
                b0count++;
                auto current_time = omp_get_wtime();
                if (b0count % 10000 == 0) {
                    std::cout << "Number of B0 assigned: " << b0count << " in " << (current_time-start_time) << " seconds." << std::endl;
                }
                if (b0count >= b0CAP)
                    break;
            }
        }

        if (omp_get_thread_num() == 0) {
            std::cout << "All tasks assigned." << std::endl;
        }

        #pragma omp taskwait
    }

    auto current_time = omp_get_wtime();
    std::cout << matrixGroup.groupName() << " centralizer count:\t" << group_centralizer_count << std::endl;
    std::cout << "Done in " << (current_time-start_time) << " seconds." << std::endl;

    b0in.close();
    b1in.close();
    centralizerFout.close();
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
    std_vector_functions::append_vector<Vector6f>(P_0, orbit_s, true);
    std_vector_functions::append_vector<Vector6f>(P_0, orbit_b, true);
    std_vector_functions::append_vector<Vector6f>(P_0, orbit_Dinvb, true);
    std_vector_functions::append_vector<Vector6f>(P_0, orbit_Dinvf, true);

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

    std_vector_functions::append_vector<Vector6f>(P_1, orbit_b, true);
    std_vector_functions::append_vector<Vector6f>(P_1, orbit_Dinvs, true);
    std_vector_functions::append_vector<Vector6f>(P_1, orbit_Dinvb, true);
    std_vector_functions::append_vector<Vector6f>(P_1, orbit_Dinvf, true);
    std_vector_functions::append_vector<Vector6f>(P_1, orbit_s, true);

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

bool AllEntriesOfParticularValues(const Matrix6f& matrix, const std::vector<float>& valid_values) {
    bool entry_is_valid;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            entry_is_valid = false;

            // check if it is any of the valid values
            for (float val : valid_values) {
                if (FloatsAreApproxEqual(matrix(i,j), val))
                    entry_is_valid = true;
            }

            if (!entry_is_valid)
                return false;
        }
    }

    // every entry is valid
    return true;
}

bool FloatsAreApproxEqual(float x, float y) {
    const float relative_difference_factor = 0.0001;    // 0.01%
    const float greater_magnitude = std::max(std::fabs(x),std::fabs(y));

    return std::fabs(x-y) < relative_difference_factor * greater_magnitude;
}