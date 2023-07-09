#include "outputResults.hpp"

int main(int argc, char *argv[]) {
    std::string curr_directory;
    // if second argument given, assume all files are in directory specified by the second argument
    if (argc > 1) {
        curr_directory = argv[1];
        if (curr_directory.back() != '/') {
            curr_directory = curr_directory + "/";
        }
    }

    std::string current_virus;
    std::string centralizer_to_check;
    std::string line;
    int b0_cols;

    std::cout << "Enter which virus to get results of:\n";
    std::getline(std::cin, current_virus);
    std::cout << "Enter which centralizer to get results of:\n";
    std::getline(std::cin, centralizer_to_check);
    std::cout << "Enter how many columns the B0 matrices will be made with: \n";
    std::getline(std::cin, line);
    b0_cols = std::stoi(line);

    if (b0_cols < 0 || b0_cols > 6) {
        std::cout << "Invalid column number input, aborting...\n";
        return 0;
    }

    const std::string OUTPUT_FILE_NAME = curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + "_" + std::to_string(b0_cols) + "_cols.txt";
    std::cout << "Getting T and B0 matrices that worked for group " + centralizer_to_check + " on virus " + current_virus + " from file...\n";
    std::ifstream fin (OUTPUT_FILE_NAME);

    if (!fin.is_open()) {
        std::cout << "Failed to open file:\t" + curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt" << std::endl;
        return 0;
    }

    outputResults::outputXcolB0(fin, b0_cols);
    fin.close();
    std::cout << "Done with file.\n";
}