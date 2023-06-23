#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <omp.h>

#include "EigenTypes.hpp"
#include "IcosahedralGroup.hpp"
#include "Matrix6fFileReader.hpp"
#include "MatrixFunctions.hpp"
#include "PermutationGroup.hpp"
#include "TetrahedralGroup.hpp"
#include "std_vector_functions.hpp"

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

using namespace EigenType;

std::vector<Vector6f> startingGeneratorsOfTCV();
std::vector<Vector6f> endingGeneratorsOfTCV();
std::vector<std::vector<Vector6f>> createOrbitsFromGenerators(const std::vector<Vector6f>& generators, Matrix6fGroup& matrix_group);

std::vector<float> possibleTransitionMatrixEntriesHardCoded();
std::vector<Matrix6f> possibleTransitionMatricesInICO(const std::vector<float>& possible_entries);
std::vector<Matrix6f> reducePossibleMatricesByEntryChecking(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& point_cloud, const std::vector<float>& valid_entries);

std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f& transition_matrix, const std::vector<std::vector<Vector6f>>& starting_orbits, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud);

void fileOutputAllMatricesOfDesiredRank(const std::vector<std::vector<Vector6f>> &P, const std::string &filename, int desiredRank, bool parallelize);

bool AllEntriesOfParticularValues(const Eigen::MatrixXf& matrix, const std::vector<float>& valid_values);
bool FloatsAreApproxEqual(float x, float y);

int main(int argc, char *argv[]) {
    std::string curr_directory;
    // if second argument given, assume all files are in directory specified by the second argument
    if (argc > 1) {
        curr_directory = argv[1];
        if (curr_directory.back() != '/') {
            curr_directory = curr_directory + "/";
        }
    }
    IcosahedralGroup icosahedral_group;

    std::vector<Vector6f> starting_generators, ending_generators;

    // eventually these should come from some sort of user input...
    starting_generators = startingGeneratorsOfTCV();
    ending_generators = endingGeneratorsOfTCV();

    std::vector<std::vector<Vector6f>> starting_orbits, ending_orbits;
    std::vector<Vector6f> starting_point_cloud, ending_point_cloud;

    Matrix6fGroup* current_matrix_group = &icosahedral_group;
    starting_orbits = createOrbitsFromGenerators(starting_generators, *current_matrix_group);
    ending_orbits = createOrbitsFromGenerators(ending_generators, *current_matrix_group);
    starting_point_cloud = std_vector_functions::unravelTwoDimVector(starting_orbits, true);
    ending_point_cloud = std_vector_functions::unravelTwoDimVector(ending_orbits, true);

    std::string current_virus = "TCV";
    std::string b0_filename = current_virus + "_B0_Matrices.csv";
    std::string b1_filename = current_virus + "_B1_Matrices.csv";


    // implement this later as it is not needed for our strategy.
    // it IS needed to figure out what the possible entries of the transition matrix can be.
    /*
    // later add ways of turning this off if desired or if the matrices have already been generated...
    int desired_matrix_rank = 6;
    FileOutputMatricesWithColumnsFromEachOrbit(starting_orbits, b0_filename);
    FileOutputMatricesWithColumnsFromEachOrbit(ending_orbits, b1_filename);
    */

    std::vector<float> possible_transition_matrix_entries = possibleTransitionMatrixEntriesHardCoded();
    std::vector<Matrix6f> possible_transition_matrices = possibleTransitionMatricesInICO(possible_transition_matrix_entries);

    std::vector<float> valid_ending_entries {-1, -0.5, 0, 0.5, 1};
    possible_transition_matrices = reducePossibleMatricesByEntryChecking(possible_transition_matrices, starting_point_cloud, valid_ending_entries);

    for (const Matrix6f& T : possible_transition_matrices) {
        std::cout << T << std::endl << std::endl;
    }
    std::cout << "Number of candidate transition matrices: " << possible_transition_matrices.size() << std::endl << std::endl;

    std::ofstream fout ("thingsThatMightWork.txt");
    std::vector<Matrix6f> possible_B0_matrices;
    Matrix6f b0_matrix;
    std::cout << "Now checking whether any potential T and B0 pairs exist..." << std::endl;
    for (const Matrix6f& transition_matrix : possible_transition_matrices) {
//        if (transition_matrix.isApprox(Matrix6f::Identity()) || transition_matrix.isApprox(Matrix6f::Identity()*-1)) {
        if (MatrixFunctions::matrixIsPlusMinusIdentity(transition_matrix)) {
            fout << transition_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
            fout << "All B0 Matrices work." << std::endl;
            fout << std::endl;
            std::cout << "Scalar multiple of identity, skipping..." << std::endl;
            continue;
        }

        possible_B0_matrices = findPossibleB0Matrices(transition_matrix, starting_orbits, starting_point_cloud, ending_point_cloud);
        std::cout << possible_B0_matrices.size() << std::endl;
        fout << transition_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
        fout << possible_B0_matrices.size() << std::endl;
        fout << std::endl;
    }

    fout.close();
}

std::vector<Vector6f> startingGeneratorsOfTCV() {
    Matrix6f D;
    D << 1,1,-1,-1,1,1,
            1,1,1,-1,-1,1,
            -1,1,1,1,-1,1,
            -1,-1,1,1,1,1,
            1,-1,-1,1,1,1,
            1,1,1,1,1,1;
    D *= 0.5;

    Vector6f s, b, f;
    s << 1, 0, 0, 0, 0, 0;
    b << 1,-1,1,1,-1,1;
    b *= 0.5;
    f << 1,0,0,-1,0,0;
    f *= 0.5;

    std::vector<Vector6f> generators;
    generators.push_back(s);
    generators.push_back(b);
    generators.emplace_back(D.inverse() * b);
    generators.emplace_back(D.inverse() * f);

    return generators;
}

std::vector<Vector6f> endingGeneratorsOfTCV() {
    Matrix6f D;
    D << 1,1,-1,-1,1,1,
            1,1,1,-1,-1,1,
            -1,1,1,1,-1,1,
            -1,-1,1,1,1,1,
            1,-1,-1,1,1,1,
            1,1,1,1,1,1;
    D *= 0.5;

    Vector6f s, b, f;
    s << 1, 0, 0, 0, 0, 0;
    b << 1,-1,1,1,-1,1;
    b *= 0.5;
    f << 1,0,0,-1,0,0;
    f *= 0.5;

    std::vector<Vector6f> generators;
    generators.push_back(s);
    generators.push_back(b);
    generators.emplace_back(D.inverse() * s);
    generators.emplace_back(D.inverse() * b);
    generators.emplace_back(D.inverse() * f);

    return generators;
}

std::vector<std::vector<Vector6f>> createOrbitsFromGenerators(const std::vector<Vector6f>& generators, Matrix6fGroup& matrix_group) {
    std::vector<std::vector<Vector6f>> orbits;
    for (const Vector6f& generator : generators) {
        orbits.push_back(matrix_group.getOrbitOfVector(generator));
    }
    return orbits;
}

std::vector<float> possibleTransitionMatrixEntriesHardCoded() {
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

std::vector<Matrix6f> possibleTransitionMatricesInICO(const std::vector<float>& possible_entries) {
    std::vector<Matrix6f> T_matrices_in_ICO;
    for (float f1 : possible_entries) {
        for (float f2: possible_entries) {
            T_matrices_in_ICO.push_back(IcosahedralGroup::matrixFormOfCentralizer(f1, f2));
        }
    }

    return T_matrices_in_ICO;
}

std::vector<Matrix6f> reducePossibleMatricesByEntryChecking(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& point_cloud, const std::vector<float>& valid_entries) {
    std::vector<Matrix6f> reduced_matrices;

    for (Matrix6f transition_matrix : possible_matrices) {
        for (const Vector6f& vector : point_cloud) {
            if (MatrixFunctions::entriesOfMatrixAreOfParticularValues(transition_matrix * vector, valid_entries)) {
                // add to list without duplicating
                MatrixFunctions::fixZeroEntries(transition_matrix);
                if (std::find(reduced_matrices.begin(), reduced_matrices.end(), transition_matrix)  == reduced_matrices.end()) {
                    reduced_matrices.push_back(transition_matrix);
                }
            }
        }
    }

    return reduced_matrices;
}

std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f& transition_matrix, const std::vector<std::vector<Vector6f>>& starting_orbits, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> possible_b0_matrices;

    std::vector<std::vector<Vector6f>> possible_columns;
    for (const std::vector<Vector6f>& orbit : starting_orbits) {
        possible_columns.push_back(orbit);
    }
    while (possible_columns.size() < 6) {
        possible_columns.push_back(starting_point_cloud);
    }

    // for each orbit, remove elements where Tm_i is not in ending point cloud
    for (std::vector<Vector6f>& column : possible_columns) {
        for (auto it = column.begin(); it != column.end();) {
            // check if we are in the ending point cloud
            bool inEndingPointCloud = false;
            for (const Vector6f& p_1 : ending_point_cloud) {
                if (p_1.isApprox(transition_matrix * (*it))) {
                    inEndingPointCloud = true;
                }
            }

            // remove element if not in ending point cloud
            if (!inEndingPointCloud) {
                it = column.erase(it);
            }
            else {
                it++;
            }
        }
    }

    // with remaining vectors in each column, brute force generate B0 matrices
    Matrix6f b0_matrix;
    for (int first = 0; first < possible_columns[0].size(); ++first) {
        for (int second = 0; second < possible_columns[1].size(); ++second) {
            for (int third = 0; third < possible_columns[2].size(); ++third) {
                for (int fourth = 0; fourth < possible_columns[3].size(); ++fourth) {
                    for (int fifth = 0; fifth < possible_columns[4].size(); ++fifth) {
                        for (int sixth = 0; sixth < possible_columns[5].size(); ++sixth) {
                            b0_matrix << possible_columns[0][first], possible_columns[1][second], possible_columns[2][third], possible_columns[3][fourth], possible_columns[4][fifth], possible_columns[5][sixth];

                            if (b0_matrix.colPivHouseholderQr().rank() == b0_matrix.cols()) {
                                possible_b0_matrices.push_back(b0_matrix);
                            }
                        }
                    }
                }
            }
        }
    }

    return possible_b0_matrices;
}