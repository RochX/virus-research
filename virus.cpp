#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <omp.h>

#include "DihedralGroupOnTen.hpp"
#include "DihedralGroupOnSix.hpp"
#include "EigenTypes.hpp"
#include "GeneratingVectorsForViruses.hpp"
#include "IcosahedralGroup.hpp"
#include "Matrix6fFileReader.hpp"
#include "MatrixFunctions.hpp"
#include "outputResults.hpp"
#include "PermutationGroup.hpp"
#include "TetrahedralGroup.hpp"
#include "std_vector_functions.hpp"

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

using namespace EigenType;

std::vector<float> possibleTransitionMatrixEntriesHardCoded(float lower_bound = -10, float upper_bound = 10);
std::vector<Matrix6f> possibleTransitionMatricesInICO(const std::vector<float>& possible_entries);
std::vector<Matrix6f> possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(const std::vector<float>& possible_entries, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud);
std::vector<Matrix6f> possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(const std::vector<float>& possible_entries, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud);
std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud);

std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                       const std::vector<Vector6f> &starting_point_cloud,
                       const std::vector<Vector6f> &ending_point_cloud);
std::vector<Matrix6fx3f> findPossibleB0Matrices3StartingOrbits(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                               const std::vector<Vector6f> &starting_point_cloud,
                                                               const std::vector<Vector6f> &ending_point_cloud);
bool verifyProductComesFromEndingOrbits(const Matrix6f& transition_matrix, const Eigen::MatrixXf& b0_matrix, const std::vector<std::vector<Vector6f>>& ending_orbits, const std::vector<Vector6f>& ending_point_cloud);

int main(int argc, char *argv[]) {
    std::string curr_directory;
    // if second argument given, assume all files are in directory specified by the second argument
    if (argc > 1) {
        curr_directory = argv[1];
        if (curr_directory.back() != '/') {
            curr_directory = curr_directory + "/";
        }
    }

    std::vector<Vector6f> starting_generators, ending_generators;
    std::string current_virus;
    std::string centralizer_to_check;
    std::string line;
    int b0_cols;

    bool get_user_input = true;
    if (get_user_input) {
        std::cout << "Enter which virus to work on:\n";
        std::getline(std::cin, current_virus);
        std::cout << "Enter which centralizer to check (ICO, A_4, D_6, D_10):\n"; // <-- this is one of: "ICO", "A_4", "D_6", or "D_10"
        std::getline(std::cin, centralizer_to_check);
        std::cout << "Enter how many columns the B0 matrices will be made with (3 or 6): \n";
        std::getline(std::cin, line);
        b0_cols = std::stoi(line);
    }
    else {
        current_virus = "1044-1127";
        centralizer_to_check = "D_10";
    }
    GeneratingVectorsForViruses::pickVirusType(current_virus, starting_generators, ending_generators);
    if (starting_generators.empty() || ending_generators.empty()) {
        std::cout << "Virus inputted is either invalid or not implemented, aborting..." << std::endl;
        return 0;
    }
    assert(!starting_generators.empty());
    assert(!ending_generators.empty());

    IcosahedralGroup icosahedral_group;
    std::vector<std::vector<Vector6f>> starting_orbits, ending_orbits;
    std::vector<Vector6f> starting_point_cloud, ending_point_cloud;

    starting_orbits = icosahedral_group.createOrbitsFromGenerators(starting_generators);
    ending_orbits = icosahedral_group.createOrbitsFromGenerators(ending_generators);
    starting_point_cloud = std_vector_functions::unravelTwoDimVector(starting_orbits, true);
    ending_point_cloud = std_vector_functions::unravelTwoDimVector(ending_orbits, true);

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

    // get all the possible ending entries
    std::vector<float> valid_ending_entries;
    for (const Vector6f& v : ending_point_cloud) {
        std::vector<float> v_entries = MatrixFunctions::entriesOfMatrix(v);
        std_vector_functions::append_vector<float>(valid_ending_entries, v_entries, true);
    }

    std::vector<float> possible_transition_matrix_entries = possibleTransitionMatrixEntriesHardCoded();
    std::vector<Matrix6f> possible_transition_matrices;

    if (centralizer_to_check == "ICO") {
        possible_transition_matrices = possibleTransitionMatricesInICO(possible_transition_matrix_entries);
        possible_transition_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(
                possible_transition_matrices,
                starting_point_cloud,
                ending_point_cloud);
    }
    else if (centralizer_to_check == "D_6") {
        possible_transition_matrices = possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(
                possible_transition_matrix_entries, starting_point_cloud, ending_point_cloud);
    }
    else if (centralizer_to_check == "D_10") {
        possible_transition_matrices = possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(
                possible_transition_matrix_entries, starting_point_cloud, ending_point_cloud);
    }
    else {
        std::cout << "Centralizer to check is of a group not implemented, aborting..." << std::endl;
        return 0;
    }

    std::cout << "First 10 candidates:" << std::endl;
    int count = 0;
    for (const Matrix6f& T : possible_transition_matrices) {
        std::cout << T << std::endl << std::endl;
        count++;
        if (count >= 10)
            break;
    }
    std::cout << "Number of candidate transition matrices: " << possible_transition_matrices.size() << std::endl << std::endl;

    std::ofstream fout (curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt");
    // for 3 column B0s
    std::vector<Matrix6fx3f> possible_B0_matrices_3_col;

    // for 6 column B0s
    std::vector<Matrix6f> possible_B0_matrices_6_col;

    std::cout << "Now checking whether any potential T and B0 pairs exist..." << std::endl;
    count = 0;
    auto start_time = omp_get_wtime();
    #pragma omp parallel private(possible_B0_matrices_3_col, possible_B0_matrices_6_col) shared(b0_cols, start_time, starting_point_cloud, ending_point_cloud, starting_orbits, ending_orbits, possible_transition_matrices, std::cout, fout, EigenType::COMMA_SEP_VALS, EigenType::TAB_INDENT) default(none)
    {
        #pragma omp master
        {
            int count = 0;
            for (const Matrix6f &transition_matrix: possible_transition_matrices) {
//        if (transition_matrix.isApprox(Matrix6f::Identity()) || transition_matrix.isApprox(Matrix6f::Identity()*-1)) {
                if (MatrixFunctions::matrixIsPlusMinusIdentity(transition_matrix)) {
                    count++;
                    std::cout << count << ":\tPlus minus identity matrix, skipping..." << std::endl;
                    continue;
                }

                // skip matrix if determinant is zero
                if (transition_matrix.determinant() < 0.000001) {
                    count++;
                    continue;
                }


                count++;
                if (b0_cols == 6) {
                    #pragma omp task private(possible_B0_matrices_6_col) firstprivate(count, transition_matrix) shared(starting_orbits, starting_point_cloud, ending_orbits, ending_point_cloud, fout, EigenType::COMMA_SEP_VALS, EigenType::TAB_INDENT, std::cout) default(none)
                    {
                        possible_B0_matrices_6_col = findPossibleB0Matrices(transition_matrix, starting_orbits,
                                                                            starting_point_cloud,
                                                                            ending_point_cloud);
                        if (!possible_B0_matrices_6_col.empty()) {
                            // check if any of the products work
                            for (const EigenType::Matrix6f &curr_b0_matrix: possible_B0_matrices_6_col) {
                                if (verifyProductComesFromEndingOrbits(transition_matrix, curr_b0_matrix, ending_orbits,
                                                                       ending_point_cloud)) {
                                    #pragma omp critical
                                    {
                                        fout << count << std::endl;
                                        fout << transition_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
                                        fout << curr_b0_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
                                        fout << (transition_matrix * curr_b0_matrix).format(EigenType::COMMA_SEP_VALS)
                                             << std::endl;
                                        fout << std::endl;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                } // endif (b0_cols == 6)
                else if (b0_cols == 3) {
                    #pragma omp task private(possible_B0_matrices_3_col) firstprivate(count, transition_matrix) shared(starting_orbits, starting_point_cloud, ending_orbits, ending_point_cloud, fout, EigenType::COMMA_SEP_VALS, EigenType::TAB_INDENT, std::cout) default(none)
                    {
                        possible_B0_matrices_3_col = findPossibleB0Matrices3StartingOrbits(transition_matrix, starting_orbits,
                                                                            starting_point_cloud,
                                                                            ending_point_cloud);
                        if (!possible_B0_matrices_3_col.empty()) {
                            // check if any of the products work
                            for (const EigenType::Matrix6fx3f &curr_b0_matrix: possible_B0_matrices_3_col) {
                                if (verifyProductComesFromEndingOrbits(transition_matrix, curr_b0_matrix, ending_orbits,
                                                                       ending_point_cloud)) {
                                    #pragma omp critical
                                    {
                                        fout << count << std::endl;
                                        fout << transition_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
                                        fout << curr_b0_matrix.format(EigenType::COMMA_SEP_VALS) << std::endl;
                                        fout << (transition_matrix * curr_b0_matrix).format(EigenType::COMMA_SEP_VALS)
                                             << std::endl;
                                        fout << std::endl;
                                    }
                                    break;
                                }
                            }
                        }
                    }
                } // endif (b0_cols == 3)

                if (count % 1000 == 0) {
                    std::cout << count << ":\tAssigned in " << omp_get_wtime()-start_time << " seconds." << std::endl;
                }
            }
        }

        #pragma omp taskwait
    }

    fout.close();

    std::cout << "Getting T and B0 matrices that worked for group " + centralizer_to_check + "...\n";
    std::ifstream fin (curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt");
    switch (b0_cols) {
        case 3:
            outputResults::output3colB0(fin);
            break;
        case 6:
            outputResults::output6colB0(fin);
            break;
        default:
            std::cout << "Invalid B0 column input." << std::endl;
    }
    std::cout << "Done with file.\n";
    fin.close();

    std::cout << "Completely done in " << omp_get_wtime()-start_time << " seconds." << std::endl;
}


std::vector<float> possibleTransitionMatrixEntriesHardCoded(float lower_bound, float upper_bound) {
    std::vector<float> possible_T_entries;

    float fourths = lower_bound;
    while (fourths <= upper_bound) {
        std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, fourths, 0.0001);
        fourths += 1.0f/4;
    }

    float sixths = lower_bound;
    while (sixths <= upper_bound) {
        std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, sixths, 0.0001);
        sixths += 1.0f/6;
    }

    std_vector_functions::push_backIfNotInVector<float>(possible_T_entries, upper_bound, 0.0001);

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

std::vector<Matrix6f> possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(const std::vector<float>& possible_entries, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> possible_transition_matrices;
    // x, y, z, t
    auto start_time = omp_get_wtime();
    std::vector<Matrix6f> reduced_matrices_first_four_vars, partial_reduced_matrices;
    int count = 0;
    Matrix6f possible;
    for (float x : possible_entries) {
        possible_transition_matrices.clear();
        for (float y: possible_entries) {
            for (float z: possible_entries) {
                for (float t: possible_entries) {
                    possible_transition_matrices.push_back(
                            DihedralGroupOnTen::matrixFormOfCentralizer(x, y, z, t, 0, 0));
                }
            }
        }

        partial_reduced_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(possible_transition_matrices,
                                                                                           starting_point_cloud,
                                                                                           ending_point_cloud);
        std_vector_functions::append_vector<Matrix6f>(reduced_matrices_first_four_vars, partial_reduced_matrices, false);

        count++;
        std::cout << "D10 checked " << count * possible_entries.size() * possible_entries.size() * possible_entries.size()
                  << " out of " << possible_entries.size() * possible_entries.size() * possible_entries.size() *
                                   possible_entries.size() << " in " << omp_get_wtime() - start_time << " seconds."
                  << std::endl;
    }

    // u, w
    std::vector<Matrix6f> reduced_matrices_last_two_vars;
    possible_transition_matrices.clear();
    for (float u : possible_entries) {
        for (float w : possible_entries) {
            possible_transition_matrices.push_back(DihedralGroupOnTen::matrixFormOfCentralizer(0,0,0,0, u, w));
        }
    }
    reduced_matrices_last_two_vars = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(
            possible_transition_matrices, starting_point_cloud, ending_point_cloud);

    // combine the two lists, hopefully the reduction did enough so this nested loop doesn't take long
    possible_transition_matrices.clear();
    for (const Matrix6f& m1 : reduced_matrices_first_four_vars) {
        for (const Matrix6f& m2 : reduced_matrices_last_two_vars) {
            possible_transition_matrices.emplace_back(m1+m2);
        }
    }

    return possible_transition_matrices;
}

std::vector<Matrix6f> possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(const std::vector<float>& possible_entries, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> possible_transition_matrices;
    // x, u, w, s
    auto start_time = omp_get_wtime();
    std::vector<Matrix6f> reduced_matrices_first_four_vars, partial_reduced_matrices;
    int count = 0;
    Matrix6f possible;
    for (float x : possible_entries) {
        possible_transition_matrices.clear();
        for (float u: possible_entries) {
            for (float w: possible_entries) {
                for (float s: possible_entries) {
                    possible_transition_matrices.push_back(
                            DihedralGroupOnSix::matrixFormOfCentralizer(x, 0, 0, 0, u, w, 0, s));
                }
            }
        }

        partial_reduced_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(possible_transition_matrices,
                                                                                           starting_point_cloud,
                                                                                           ending_point_cloud);

        std_vector_functions::append_vector<Matrix6f>(reduced_matrices_first_four_vars, partial_reduced_matrices, false);

        count++;
        std::cout << "First half D6 checked " << count * possible_entries.size() * possible_entries.size() * possible_entries.size()
                  << " out of " << possible_entries.size() * possible_entries.size() * possible_entries.size() *
                                   possible_entries.size() << " in " << omp_get_wtime() - start_time << " seconds."
                  << std::endl;
    }

    // y, z, t, v
    count = 0;
    std::vector<Matrix6f> reduced_matrices_last_four_vars;
    possible_transition_matrices.clear();
    for (float y : possible_entries) {
        possible_transition_matrices.clear();
        for (float z: possible_entries) {
            for (float t: possible_entries) {
                for (float v: possible_entries) {
                    possible_transition_matrices.push_back(
                            DihedralGroupOnSix::matrixFormOfCentralizer(0, y, z, t, 0, 0, v, 0));
                }
            }
        }

        partial_reduced_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(possible_transition_matrices,
                                                                                           starting_point_cloud,
                                                                                           ending_point_cloud);

        std_vector_functions::append_vector<Matrix6f>(reduced_matrices_last_four_vars, partial_reduced_matrices, false);

        count++;
        std::cout << "Second half D6 checked " << count * possible_entries.size() * possible_entries.size() * possible_entries.size()
                  << " out of " << possible_entries.size() * possible_entries.size() * possible_entries.size() *
                                   possible_entries.size() << " in " << omp_get_wtime() - start_time << " seconds."
                  << std::endl;
    }

    // combine the two lists, hopefully the reduction did enough so this nested loop doesn't take long
    possible_transition_matrices.clear();
    for (const Matrix6f& m1 : reduced_matrices_first_four_vars) {
        for (const Matrix6f& m2 : reduced_matrices_last_four_vars) {
            possible_transition_matrices.emplace_back(m1+m2);
        }
    }

    return possible_transition_matrices;
}

std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> reduced_matrices;
    float max_norm_of_ending_point_cloud = 0;

    // find max norm of ending point cloud
    for (const Vector6f& v : ending_point_cloud) {
        max_norm_of_ending_point_cloud = std::max(max_norm_of_ending_point_cloud, v.norm());
    }

    int count = 0;
    auto start_time = omp_get_wtime();
    Vector6f current_Tv_product;
    for (Matrix6f transition_matrix : possible_matrices) {
        for (const Vector6f& starting_vector : starting_point_cloud) {
            current_Tv_product = transition_matrix*starting_vector;
            // we can skip this pair if T*v has norm greater than any element of the ending point cloud.
            if (current_Tv_product.norm() > max_norm_of_ending_point_cloud)
                continue;

            // add to list, no need to check for duplicates because we get each unique tuple of variables.
            MatrixFunctions::fixZeroEntries(transition_matrix);
            reduced_matrices.push_back(transition_matrix);
            break;
        }
    }

    return reduced_matrices;
}

std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                       const std::vector<Vector6f> &starting_point_cloud,
                       const std::vector<Vector6f> &ending_point_cloud) {
    std::vector<Matrix6f> possible_b0_matrices;

    std::vector<std::vector<Vector6f>> possible_columns;
    for (const std::vector<Vector6f>& orbit : starting_orbits) {
        possible_columns.push_back(orbit);
    }
    while (possible_columns.size() < 6) {
        possible_columns.push_back(starting_point_cloud);
    }
    int total = 1;
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

        total *= column.size();
    }

//    if (total > 100000) {
//        std::cout << "Matrix has " << total << " which is a lot to check..." << std::endl;
//        std::cout << transition_matrix << std::endl;
//    }

    // with remaining vectors in each column, brute force generate B0 matrices
    Matrix6f b0_matrix;
    b0_matrix.setZero();
    for (int first = 0; first < possible_columns[0].size(); ++first) {
        b0_matrix.col(0) << possible_columns[0][first];

        for (int second = 0; second < possible_columns[1].size(); ++second) {
            b0_matrix.col(1) << possible_columns[1][second];
            if (b0_matrix.colPivHouseholderQr().rank() < 2)
                continue;

            for (int third = 0; third < possible_columns[2].size(); ++third) {
                b0_matrix.col(2) << possible_columns[2][third];
                if (b0_matrix.colPivHouseholderQr().rank() < 3)
                    continue;

                for (int fourth = 0; fourth < possible_columns[3].size(); ++fourth) {
                    b0_matrix.col(3) << possible_columns[3][fourth];
                    if (b0_matrix.colPivHouseholderQr().rank() < 4)
                        continue;

                    for (int fifth = 0; fifth < possible_columns[4].size(); ++fifth) {
                        b0_matrix.col(4) << possible_columns[4][fifth];
                        if (b0_matrix.colPivHouseholderQr().rank() < 5)
                            continue;

                        for (int sixth = 0; sixth < possible_columns[5].size(); ++sixth) {
                            b0_matrix.col(5) << possible_columns[5][sixth];

                            if (b0_matrix.colPivHouseholderQr().rank() == b0_matrix.cols()) {
                                possible_b0_matrices.push_back(b0_matrix);
                            }

                            b0_matrix.col(5).setZero();
                        }
                        b0_matrix.col(4).setZero();
                    }
                    b0_matrix.col(3).setZero();
                }
                b0_matrix.col(2).setZero();
            }
            b0_matrix.col(1).setZero();
        }
    }

    return possible_b0_matrices;
}

std::vector<Matrix6fx3f> findPossibleB0Matrices3StartingOrbits(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                                             const std::vector<Vector6f> &starting_point_cloud,
                                             const std::vector<Vector6f> &ending_point_cloud) {
    assert(starting_orbits.size() == 3);

    std::vector<Matrix6fx3f> possible_b0_matrices;

    std::vector<std::vector<Vector6f>> possible_columns;
    for (const std::vector<Vector6f>& orbit : starting_orbits) {
        possible_columns.push_back(orbit);
    }

    unsigned long long total = 1;
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

        total *= column.size();
    }

    Matrix6fx3f b0_matrix;
    b0_matrix.setZero();
    for (int first = 0; first < possible_columns[0].size(); ++first) {
        for (int second = 0; second < possible_columns[1].size(); ++second) {
            for (int third = 0; third < possible_columns[2].size(); ++third) {
                b0_matrix.col(0) = possible_columns[0][first];
                b0_matrix.col(1) = possible_columns[1][second];
                b0_matrix.col(2) = possible_columns[2][third];

                if (b0_matrix.colPivHouseholderQr().rank() == b0_matrix.cols())
                    possible_b0_matrices.push_back(b0_matrix);
            }
        }
    }

    return possible_b0_matrices;
}


bool verifyProductComesFromEndingOrbits(const Matrix6f& transition_matrix, const Eigen::MatrixXf& b0_matrix, const std::vector<std::vector<Vector6f>>& ending_orbits, const std::vector<Vector6f>& ending_point_cloud) {
    Eigen::MatrixXf product = transition_matrix*b0_matrix;
    // check for full rank
    if (product.colPivHouseholderQr().rank() < product.cols())
        return false;

    // check each column actually comes from the ending point cloud
    auto it = ending_point_cloud.begin();
    for (int i = 0; i < product.cols(); i++) {
        for (it = ending_point_cloud.begin(); it != ending_point_cloud.end(); it++) {
            if (product.col(i).isApprox(*it))
                break;
        }

        // in current column we are not from the ending point cloud...
        if (it == ending_point_cloud.end())
            return false;
    }

    // check for a representative from each orbit
    int orbits_represented = 0;
    bool curr_orbit_done;
    for (const std::vector<Vector6f>& orbit : ending_orbits) {
        curr_orbit_done = false;
        for (const Vector6f& vector : orbit) {
            for (int i = 0; i < product.cols(); i++) {
                if (vector.isApprox(product.col(i))) {
                    orbits_represented++;
                    curr_orbit_done = true;
                    break;
                }
            }

            if (curr_orbit_done)
                break;
        }
    }

    return orbits_represented == ending_orbits.size();
}