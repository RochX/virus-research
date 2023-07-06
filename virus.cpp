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
std::vector<Matrix6f> possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                   const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                   const std::vector<Vector6f> &starting_point_cloud,
                                                                   const std::vector<Vector6f> &ending_point_cloud);
std::vector<Matrix6f> possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                  const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                  const std::vector<Vector6f> &starting_point_cloud,
                                                                  const std::vector<Vector6f> &ending_point_cloud);
std::vector<Matrix6f> possibleTransitionMatricesInA4WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                  const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                  const std::vector<Vector6f> &starting_point_cloud,
                                                                  const std::vector<Vector6f> &ending_point_cloud);
// TODO: add user input to choose which one of these to use (i.e. full check or partial check)
// for now use the full check
std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud);
std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<std::vector<Vector6f>>& starting_orbits, const std::vector<Vector6f>& ending_point_cloud);

std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                       const std::vector<Vector6f> &starting_point_cloud,
                       const std::vector<Vector6f> &ending_point_cloud);
std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f &transition_matrix, const std::vector<std::vector<Vector6f>> &starting_orbits,
                                             const std::vector<Vector6f> &starting_point_cloud,
                                             const std::vector<Vector6f> &ending_point_cloud, int num_cols);
bool verifyProductComesFromEndingOrbits(const Matrix6f& transition_matrix, const Eigen::MatrixXf& b0_matrix, const std::vector<std::vector<Vector6f>>& ending_orbits, const std::vector<Vector6f>& ending_point_cloud, int num_cols = -1);

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

    std::cout << "Enter which virus to work on:\n";
    std::getline(std::cin, current_virus);
    GeneratingVectorsForViruses::pickVirusType(current_virus, starting_generators, ending_generators);
    std::cout
            << "Enter which centralizer to check (ICO, A_4, D_6, D_10):\n";// <-- this is one of: "ICO", "A_4", "D_6", or "D_10"
    GeneratingVectorsForViruses::pickVirusType(current_virus, starting_generators, ending_generators);
    std::getline(std::cin, centralizer_to_check);
    std::cout << "The current virus has " << starting_generators.size() << " starting generators and "
              << ending_generators.size() << " ending generators.\n";
    std::cout << "Enter how many columns the B0 matrices will be made with:\n";
    std::cout << "(for good results, use something between " << std::max(starting_generators.size(), ending_generators.size()) <<" and 6, inclusive)\n";
    std::getline(std::cin, line);
    b0_cols = std::stoi(line);
    if (starting_generators.size() != ending_generators.size()) {
        std::cout << "Unequal number of starting and ending generators, swap the start and end generators?\n";
        std::cout << "Swapping when there are more ending generators than starting generators usually leads to faster runtime.\n";
        std::cout << "(Y/N):\t";
        std::getline(std::cin, line);
        if (line == "Y" || line == "y") {
            std::cout << "Swapping, note that you will now need to take the inverse of any transition this program finds."
                      << std::endl;
            starting_generators.swap(ending_generators);
        }
    }
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

    // TODO: add more versatile user input, for now this is fine
    std::vector<float> possible_transition_matrix_entries = possibleTransitionMatrixEntriesHardCoded(-5, 5);
    std::vector<std::vector<float>> example_entry_lists {{-1, 0, 1},
                                                         {-1, -0.5, 0, 0.5, 1},
                                                         {-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2},
                                                         {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3},
                                                         possibleTransitionMatrixEntriesHardCoded(-1, 1),
                                                         possibleTransitionMatrixEntriesHardCoded(-2, 2),
                                                         possibleTransitionMatrixEntriesHardCoded(-3, 3)};
    std::cout << "Which entry list for transition matrices to use? Here are the lists:\n";
    for (int i = 0; i < example_entry_lists.size(); i++) {
        std::cout << "List " << i << ":\t";
        for (float f : example_entry_lists[i]) {
            std::cout << f;
            if (f != example_entry_lists[i].back())
                std::cout << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "Now pick one of the above lists:\n";
    getline(std::cin, line);
    int user_list_pick = std::stoi(line);
    if (0 <= user_list_pick && user_list_pick < example_entry_lists.size()) {
        possible_transition_matrix_entries = example_entry_lists[user_list_pick];
    }
    else {
        std::cout << "Invalid list picked, aborting...\n";
        return 0;
    }

    auto program_start_time = omp_get_wtime();
    std::vector<Matrix6f> possible_transition_matrices;

    if (centralizer_to_check == "ICO") {
        possible_transition_matrices = possibleTransitionMatricesInICO(possible_transition_matrix_entries);
        possible_transition_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(
                possible_transition_matrices,
                starting_point_cloud,
                ending_point_cloud);
    }
    else if (centralizer_to_check == "A_4") {
        possible_transition_matrices = possibleTransitionMatricesInA4WithCheckingMapIntoEndingPointCloud(
                possible_transition_matrix_entries, starting_orbits, starting_point_cloud, ending_point_cloud);
    }
    else if (centralizer_to_check == "D_6") {
        possible_transition_matrices = possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(
                possible_transition_matrix_entries, starting_orbits, starting_point_cloud, ending_point_cloud);
    }
    else if (centralizer_to_check == "D_10") {
        possible_transition_matrices = possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(
                possible_transition_matrix_entries, starting_orbits, starting_point_cloud, ending_point_cloud);
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

    std::vector<Matrix6f> possible_B0_matrices;

    std::cout << "Now checking whether any potential T and B0 pairs exist..." << std::endl;
    count = 0;
    auto start_time = omp_get_wtime();
    #pragma omp parallel private(possible_B0_matrices) shared(b0_cols, start_time, starting_point_cloud, ending_point_cloud, starting_orbits, ending_orbits, possible_transition_matrices, std::cout, fout, EigenType::COMMA_SEP_VALS, EigenType::TAB_INDENT) default(none)
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

                #pragma omp task private(possible_B0_matrices) firstprivate(count, transition_matrix) shared(b0_cols, starting_orbits, starting_point_cloud, ending_orbits, ending_point_cloud, fout, EigenType::COMMA_SEP_VALS, EigenType::TAB_INDENT, std::cout) default(none)
                {
                    possible_B0_matrices = findPossibleB0Matrices(transition_matrix, starting_orbits,
                                                                        starting_point_cloud,
                                                                        ending_point_cloud, b0_cols);
                    if (!possible_B0_matrices.empty()) {
                        // check if any of the products work
                        for (const EigenType::Matrix6f &curr_b0_matrix: possible_B0_matrices) {
                            if (verifyProductComesFromEndingOrbits(transition_matrix, curr_b0_matrix, ending_orbits,
                                                                   ending_point_cloud, b0_cols)) {
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

                if (count % 1000 == 0) {
                    std::cout << count << ":\tAssigned in " << omp_get_wtime()-start_time << " seconds." << std::endl;
                }
            }
        }

        #pragma omp taskwait
    }

    fout.close();

    std::cout << "Getting T and B0 matrices that worked for group " + centralizer_to_check + " from file...\n";
    std::ifstream fin (curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt");
    outputResults::outputXcolB0(fin, b0_cols);
    fin.close();
    std::cout << "Done with file.\n";

    std::cout << "Completely done in " << omp_get_wtime()-program_start_time << " seconds." << std::endl;
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

std::vector<Matrix6f> possibleTransitionMatricesInA4WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                  const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                  const std::vector<Vector6f> &starting_point_cloud,
                                                                  const std::vector<Vector6f> &ending_point_cloud) {
    std::vector<Matrix6f> possible_transition_matrices;
    // x, y, z, t
    auto start_time = omp_get_wtime();
    std::vector<Matrix6f> partial_reduced_matrices;
    int count = 0;
    Matrix6f possible;
    for (float x : possible_entries) {
        for (float y: possible_entries) {
            for (float z: possible_entries) {
                for (float t: possible_entries) {
                    possible_transition_matrices.push_back(
                            TetrahedralGroup::matrixFormOfCentralizer(x, y, z, t));
                }
            }
        }

        partial_reduced_matrices = reducePossibleMatricesByCheckingMapIntoEndingPointCloud(possible_transition_matrices,
                                                                                           starting_orbits,
                                                                                           ending_point_cloud);
        std_vector_functions::append_vector<Matrix6f>(possible_transition_matrices, partial_reduced_matrices, true);

        count++;
        std::cout << "A4 checked " << count * possible_entries.size() * possible_entries.size() * possible_entries.size()
                  << " out of " << possible_entries.size() * possible_entries.size() * possible_entries.size() *
                                   possible_entries.size() << " in " << omp_get_wtime() - start_time << " seconds."
                  << std::endl;
    }

    return possible_transition_matrices;
}

std::vector<Matrix6f> possibleTransitionMatricesInD10WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                   const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                   const std::vector<Vector6f> &starting_point_cloud,
                                                                   const std::vector<Vector6f> &ending_point_cloud) {
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
                                                                                           starting_orbits,
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
            possible_transition_matrices, starting_orbits, ending_point_cloud);

    // combine the two lists, hopefully the reduction did enough so this nested loop doesn't take long
    std::cout << "Merging D_10 lists together, a total of " << reduced_matrices_first_four_vars.size()*reduced_matrices_last_two_vars.size() << " to merge..." << std::endl;
    possible_transition_matrices.clear();
    for (const Matrix6f& m1 : reduced_matrices_first_four_vars) {
        for (const Matrix6f& m2 : reduced_matrices_last_two_vars) {
            // only consider matrices that have nonzero determinant
            if (std::fabs((m1+m2).determinant()) >= 0.00001)
                possible_transition_matrices.emplace_back(m1+m2);
        }
    }
    std::cout << "Done with merging." << std::endl;

    return possible_transition_matrices;
}

std::vector<Matrix6f> possibleTransitionMatricesInD6WithCheckingMapIntoEndingPointCloud(const std::vector<float> &possible_entries,
                                                                  const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                  const std::vector<Vector6f> &starting_point_cloud,
                                                                  const std::vector<Vector6f> &ending_point_cloud) {
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
                                                                                           starting_orbits,
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
                                                                                           starting_orbits,
                                                                                           ending_point_cloud);

        std_vector_functions::append_vector<Matrix6f>(reduced_matrices_last_four_vars, partial_reduced_matrices, false);

        count++;
        std::cout << "Second half D6 checked " << count * possible_entries.size() * possible_entries.size() * possible_entries.size()
                  << " out of " << possible_entries.size() * possible_entries.size() * possible_entries.size() *
                                   possible_entries.size() << " in " << omp_get_wtime() - start_time << " seconds."
                  << std::endl;
    }

    // combine the two lists, hopefully the reduction did enough so this nested loop doesn't take long
    std::cout << "Merging D_6 lists together, a total of " << reduced_matrices_first_four_vars.size()*reduced_matrices_last_four_vars.size() << " to merge..." << std::endl;
    possible_transition_matrices.clear();
    for (const Matrix6f& m1 : reduced_matrices_first_four_vars) {
        for (const Matrix6f& m2 : reduced_matrices_last_four_vars) {
            // only consider matrices that have nonzero determinant
            if (std::fabs((m1+m2).determinant()) >= 0.00001)
                possible_transition_matrices.emplace_back(m1+m2);
        }
    }
    std::cout << "Done with merging." << std::endl;

    return possible_transition_matrices;
}

std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<Vector6f>& starting_point_cloud, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> reduced_matrices;
    float max_norm_of_ending_point_cloud = 0;
    std::vector<float> ending_entries;

    // find max norm of ending point cloud and ending entries
    for (const Vector6f& v : ending_point_cloud) {
        max_norm_of_ending_point_cloud = std::max(max_norm_of_ending_point_cloud, v.norm());
        std::vector<float> curr_entries = MatrixFunctions::entriesOfMatrix(v);
        std_vector_functions::append_vector(ending_entries, curr_entries, true);
    }
    std_vector_functions::push_backIfNotInVector<float>(ending_entries, 0, 0.0001);

    int count = 0;
    auto start_time = omp_get_wtime();
    Vector6f current_Tv_product;
    for (Matrix6f transition_matrix : possible_matrices) {
        for (const Vector6f& starting_vector : starting_point_cloud) {
            current_Tv_product = transition_matrix*starting_vector;
            // we can skip this pair if T*v has norm greater than any element of the ending point cloud.
            if (current_Tv_product.norm() > max_norm_of_ending_point_cloud)
                continue;

            // check if the entries are correct
            if (MatrixFunctions::entriesOfMatrixAreOfParticularValues(current_Tv_product, ending_entries)) {
                // add to list, no need to check for duplicates because we get each unique tuple of variables.
                MatrixFunctions::fixZeroEntries(transition_matrix);
                reduced_matrices.push_back(transition_matrix);
                break;
            }
        }
    }

    return reduced_matrices;
}

std::vector<Matrix6f> reducePossibleMatricesByCheckingMapIntoEndingPointCloud(const std::vector<Matrix6f>& possible_matrices, const std::vector<std::vector<Vector6f>>& starting_orbits, const std::vector<Vector6f>& ending_point_cloud) {
    std::vector<Matrix6f> reduced_matrices;
    float max_norm_of_ending_point_cloud = 0;

    // find max norm of ending point cloud
    for (const Vector6f& v : ending_point_cloud) {
        max_norm_of_ending_point_cloud = std::max(max_norm_of_ending_point_cloud, v.norm());
    }

    int orbits_that_might_work;
    auto start_time = omp_get_wtime();
    Vector6f current_Tv_product;
    for (Matrix6f transition_matrix : possible_matrices) {
        orbits_that_might_work = 0;
        for (const std::vector<Vector6f>& orbit : starting_orbits) {
            bool orbits_maps_into_ending_point_cloud = [&] {
                // start lambda function
                for (const Vector6f &starting_vector: orbit) {
                    // eliminate based on norm...
                    if ((transition_matrix * starting_vector).norm() > max_norm_of_ending_point_cloud) {
                        continue;
                    }

                    for (Vector6f ending_vector: ending_point_cloud) {
                        // if transition matrix has a zero row, zero out that entry of this particular vector
                        for (int i = 0; i < transition_matrix.rows(); i++) {
                            // row is (approx) zero row
                            if (transition_matrix.row(i).norm() < 0.00001) {
                                ending_vector[i] = 0;
                            }
                        }

                        if ((transition_matrix * starting_vector).isApprox(ending_vector)) {
                            return true;
                        }
                    }
                }
                return false;
            }(); // end lambda function
            if (orbits_maps_into_ending_point_cloud)
                orbits_that_might_work++;
        }
        // check if something from each orbit could map into the ending point cloud
        if (orbits_that_might_work == starting_orbits.size()) {
            MatrixFunctions::fixZeroEntries(transition_matrix);
            reduced_matrices.push_back(transition_matrix);
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

bool verifyProductComesFromEndingOrbits(const Matrix6f& transition_matrix, const Eigen::MatrixXf& b0_matrix, const std::vector<std::vector<Vector6f>>& ending_orbits, const std::vector<Vector6f>& ending_point_cloud, int num_cols) {
    Eigen::MatrixXf product = transition_matrix*b0_matrix;
    if (num_cols < 0) {
        num_cols = product.cols();
    }
    assert(0 <= num_cols);
    assert(num_cols <= 6);

    // check for full rank
    if (product.colPivHouseholderQr().rank() < num_cols)
        return false;

    // check each nonzero column actually comes from the ending point cloud
    auto it = ending_point_cloud.begin();
    for (int i = 0; i < num_cols; i++) {
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

// use a variable number of columns by just filling in the ending columns with zeroes
std::vector<Matrix6f> findPossibleB0Matrices(const Matrix6f &transition_matrix,
                                                                const std::vector<std::vector<Vector6f>> &starting_orbits,
                                                                const std::vector<Vector6f> &starting_point_cloud,
                                                                const std::vector<Vector6f> &ending_point_cloud,
                                                                int num_cols) {
    assert(num_cols <= 6);
    assert(starting_orbits.size() <= num_cols);
    std::vector<Matrix6f> possible_b0_matrices;

    std::vector<std::vector<Vector6f>> possible_columns;
    for (const std::vector<Vector6f>& orbit : starting_orbits) {
        possible_columns.push_back(orbit);
    }
    while (possible_columns.size() < num_cols) {
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

    // none to check ==> nothing to do
    if (total == 0)
        return possible_b0_matrices;

    // to make matrices 6x6, fill in remaining columns with zeroes
    std::vector<Vector6f> zero_column {Vector6f::Zero()};
    while (possible_columns.size() < 6) {
        possible_columns.push_back(zero_column);
    }

    // with remaining vectors in each column, brute force generate B0 matrices
    Matrix6f b0_matrix;
    b0_matrix.setZero();
    for (int first = 0; first < possible_columns[0].size(); ++first) {
        b0_matrix.col(0) << possible_columns[0][first];

        for (int second = 0; second < possible_columns[1].size(); ++second) {
            b0_matrix.col(1) << possible_columns[1][second];
            if (b0_matrix.colPivHouseholderQr().rank() < std::min(2, num_cols))
                continue;

            for (int third = 0; third < possible_columns[2].size(); ++third) {
                b0_matrix.col(2) << possible_columns[2][third];
                if (b0_matrix.colPivHouseholderQr().rank() < std::min(3, num_cols))
                    continue;

                for (int fourth = 0; fourth < possible_columns[3].size(); ++fourth) {
                    b0_matrix.col(3) << possible_columns[3][fourth];
                    if (b0_matrix.colPivHouseholderQr().rank() < std::min(4, num_cols))
                        continue;

                    for (int fifth = 0; fifth < possible_columns[4].size(); ++fifth) {
                        b0_matrix.col(4) << possible_columns[4][fifth];
                        if (b0_matrix.colPivHouseholderQr().rank() < std::min(5, num_cols))
                            continue;

                        for (int sixth = 0; sixth < possible_columns[5].size(); ++sixth) {
                            b0_matrix.col(5) << possible_columns[5][sixth];

                            if (b0_matrix.colPivHouseholderQr().rank() == num_cols) {
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