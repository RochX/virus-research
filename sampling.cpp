#include "EigenTypes.hpp"
#include "GeneratingVectorsForViruses.hpp"
#include "IcosahedralGroup.hpp"
#include "Matrix6fFileReader.hpp"
#include "Matrix6fGroup.hpp"
#include "MatrixFunctions.hpp"
#include "std_vector_functions.hpp"

using namespace EigenType;

std::vector<float> findPossibleTEntriesBySampling(std::ifstream &b0_in, std::ifstream &b1_in, int num_samples);
void fileOutputAllMatricesOfFullRank(const std::vector<std::vector<Vector6f>> &orbits,
                                     const std::vector<Vector6f> &point_cloud, const std::string &filename,
                                     int sample_size, bool randomize);

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
    std::vector<std::vector<Vector6f>> starting_orbits, ending_orbits;
    std::vector<Vector6f> starting_point_cloud, ending_point_cloud;
    std::string current_virus, line;

    std::cout << "Enter which virus to work on:\n";
    std::getline(std::cin, current_virus);
    GeneratingVectorsForViruses::pickVirusType(current_virus, starting_generators, ending_generators);

    if (starting_generators.empty() || ending_generators.empty()) {
        std::cout << "Virus inputted is either invalid or not implemented, aborting..." << std::endl;
        return 0;
    }
    assert(!starting_generators.empty());
    assert(!ending_generators.empty());

    std::cout << "The current virus has " << starting_generators.size() << " starting generators and "
              << ending_generators.size() << " ending generators.\n";
    if (starting_generators.size() != ending_generators.size()) {
        std::cout << "Unequal number of starting and ending generators, swap the start and end generators?\n";
        std::cout << "Swapping when there are more ending generators than starting generators usually leads to faster runtime.\n";
        std::cout << "(Y/N):";
        std::getline(std::cin, line);
        if (line == "Y" || line == "y") {
            starting_generators.swap(ending_generators);
        }
    }

    std::cout << "Input sample size:";
    std::getline(std::cin, line);
    int sample_size = std::stoi(line);
    assert(sample_size > 0);
    std::cout << "Number to divide sample size into:";
    std::getline(std::cin, line);
    int divides = std::stoi(line);


    IcosahedralGroup icosahedral_group;
    starting_orbits = icosahedral_group.createOrbitsFromGenerators(starting_generators);
    ending_orbits = icosahedral_group.createOrbitsFromGenerators(ending_generators);

    starting_point_cloud = std_vector_functions::unravelTwoDimVector(starting_orbits);
    ending_point_cloud = std_vector_functions::unravelTwoDimVector(ending_orbits);

    std::string b0_filename = curr_directory + current_virus + "_B0_matrices.csv";
    std::string b1_filename = curr_directory + current_virus + "_B1_matrices.csv";

    fileOutputAllMatricesOfFullRank(starting_orbits, starting_point_cloud, b0_filename, sample_size, true);
    fileOutputAllMatricesOfFullRank(ending_orbits, ending_point_cloud, b1_filename, sample_size, true);

    std::ifstream b0_in(b0_filename);
    std::ifstream b1_in(b1_filename);
    // read csv headers
    std::getline(b0_in, line);
    std::getline(b1_in, line);
    std::vector<float> possible_entries;

    for (int i = 0; i < sample_size/divides; i++) {
        possible_entries = findPossibleTEntriesBySampling(b0_in, b1_in, divides);
        std::cout << "Entries:\t";
        for (float e: possible_entries) {
            std::cout << e;
            if (e != possible_entries.back())
                std::cout << ", ";
        }
        std::cout << std::endl;
    }
}

std::vector<float> getTEntries(const Matrix6f& b0_matrix, const Matrix6f& b1_matrix) {
    Matrix6f T = b1_matrix * b0_matrix.inverse();
    MatrixFunctions::fixZeroEntries(T);
    return MatrixFunctions::entriesOfMatrix(T);
}

std::vector<float> findPossibleTEntriesBySampling(std::ifstream &b0_in, std::ifstream &b1_in, int num_samples) {
    std::vector<float> possible_T_entries, curr_entries;
    std::vector<Matrix6f> b0_matrices, b1_matrices;

    Matrix6fFileReader::readNextNMatrices(b0_in, b0_matrices, num_samples);
    Matrix6fFileReader::readNextNMatrices(b1_in, b1_matrices, num_samples);

    for (const Matrix6f& b0 : b0_matrices) {
        for (const Matrix6f& b1 : b1_matrices) {
            curr_entries = getTEntries(b0, b1);
            std_vector_functions::append_vector(possible_T_entries, curr_entries, true);
        }
    }

    std::sort(possible_T_entries.begin(), possible_T_entries.end());
    return possible_T_entries;
}

void fileOutputAllMatricesOfFullRank(const std::vector<std::vector<Vector6f>> &orbits,
                                     const std::vector<Vector6f> &point_cloud, const std::string &filename,
                                     int sample_size, bool randomize) {
    std::string CSV_HEADER = "11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 61, 62, 63, 64, 65, 66";


    std::vector<std::vector<Vector6f>> column_choices;
    for (const std::vector<Vector6f>& orbit : orbits) {
        column_choices.push_back(orbit);
    }
    while (column_choices.size() < 6) {
        column_choices.push_back(point_cloud);
    }

    Matrix6f m;

    // check if we opened file properly, if so output csv header
    std::ofstream fout (filename);
    if (!fout.is_open()) {
        std::cerr << "Failed to open file " << filename << " within function fileOutputAllMatricesOfFullRank." << std::endl;
        return;
    }
    else {
        fout << CSV_HEADER << std::endl;
    }

    std::cout << "Starting full rank matrix output to file:\t" + filename << std::endl;

    // start tracking time
    auto start_time = omp_get_wtime();
    long long total_full_rank_matrices = 0;

    if (randomize) {
        int random_indices[6];
        while (total_full_rank_matrices < sample_size) {
            for (int i = 0; i < 6; i++) {
                random_indices[i] = std::rand() % column_choices[i].size();
            }
            for (int i = 0; i < 6; i++) {
                m.col(i) << column_choices[i][random_indices[i]];
            }
            if (m.colPivHouseholderQr().rank() == m.cols()) {
                fout << m.format(EigenType::COMMA_SEP_VALS) << std::endl;
                total_full_rank_matrices++;
            }
        }
    }
    else {
        for (int first = 0; first < column_choices[0].size(); ++first) {
            for (int second = 0; second < column_choices[1].size(); ++second) {
                for (int third = 0; third < column_choices[2].size(); ++third) {
                    for (int fourth = 0; fourth < column_choices[3].size(); ++fourth) {
                        for (int fifth = 0; fifth < column_choices[4].size(); ++fifth) {
                            for (int sixth = 0; sixth < column_choices[5].size(); ++sixth) {
                                m
                                        << column_choices[0][first], column_choices[1][second], column_choices[2][third], column_choices[3][fourth], column_choices[4][fifth], column_choices[5][sixth];

                                // check our matrix is full rank
                                if (m.colPivHouseholderQr().rank() == m.cols()) {
                                    fout << m.format(EigenType::COMMA_SEP_VALS) << std::endl;
                                    total_full_rank_matrices++;
                                }

                                if (total_full_rank_matrices >= sample_size)
                                    goto finishOutput;
                            }
                        }
                    }
                }
            }
        }
    }
    finishOutput:
    // finished with parallel region means we finished with file
    fout.close();

    // finish tracking time, output some results
    auto current_time = omp_get_wtime();
    std::cout << "\nFinished outputting to " << filename << "." << std::endl;
    std::cout << "Outputted " << total_full_rank_matrices << " matrices of desired rank 6." << std::endl;
    std::cout << "Total time taken:\t" << (current_time - start_time) << " seconds" << std::endl << std::endl;
}

/*
 * else if (centralizer_to_check == "testing") {
        float x, y, z, t, u, w, v, s;
        std::vector<float> inverse_entries;
        Matrix6f T_inverse;
        auto print_inverse_entries_sorted = [&] {
            std::sort(inverse_entries.begin(), inverse_entries.end());
            for (float f : inverse_entries) {
                std::cout << f << " ";
            }
        };
        auto lambda = [&] {
            int random_pos = std::rand() % possible_transition_matrix_entries.size();
            return possible_transition_matrix_entries[random_pos];
        };
        std::cout << "ICO inverse entries:" << std::endl;
        for (int k = 0; k < 100; k++) {
            inverse_entries.clear();
            x = lambda();
            y = lambda();
            z = lambda();
            t = lambda();
            u = lambda();
            w = lambda();
            v = lambda();
            s = lambda();
            T_inverse = IcosahedralGroup::matrixFormOfCentralizer(x, y).inverse();
            inverse_entries = MatrixFunctions::entriesOfMatrix(T_inverse, true);
            print_inverse_entries_sorted();
        }

        std::cout << "A_4 inverse entries:" << std::endl;
        for (int k = 0; k < 100; k++) {
            inverse_entries.clear();
            x = lambda();
            y = lambda();
            z = lambda();
            t = lambda();
            u = lambda();
            w = lambda();
            v = lambda();
            s = lambda();
            T_inverse = TetrahedralGroup::matrixFormOfCentralizer(x, y, z, t).inverse();
            inverse_entries = MatrixFunctions::entriesOfMatrix(T_inverse);
            print_inverse_entries_sorted();
        }

        std::cout << "D_6 inverse entries:" << std::endl;
        for (int k = 0; k < 100; k++) {
            inverse_entries.clear();
            x = lambda();
            y = lambda();
            z = lambda();
            t = lambda();
            u = lambda();
            w = lambda();
            v = lambda();
            s = lambda();
            T_inverse = DihedralGroupOnSix::matrixFormOfCentralizer(x, y, z, t, u, w, v, s).inverse();
            inverse_entries = MatrixFunctions::entriesOfMatrix(T_inverse);
            print_inverse_entries_sorted();
        }

        std::cout << "D_10 inverse entries:" << std::endl;
        for (int k = 0; k < 100; k++) {
            inverse_entries.clear();
            x = lambda();
            y = lambda();
            z = lambda();
            t = lambda();
            u = lambda();
            w = lambda();
            v = lambda();
            s = lambda();
            T_inverse = DihedralGroupOnTen::matrixFormOfCentralizer(x, y, z, t, u, w).inverse();
            inverse_entries = MatrixFunctions::entriesOfMatrix(T_inverse);
            print_inverse_entries_sorted();
        }
 */