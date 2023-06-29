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

    bool get_user_input = true;
    if (get_user_input) {
        std::cout << "Enter which virus to get results of:\n";
        std::getline(std::cin, current_virus);
        std::cout << "Enter which centralizer to get results of:\n";
        std::getline(std::cin, centralizer_to_check);
        std::cout << "Enter how many columns the B0 matrices will be made with (3 or 6): \n";
        std::getline(std::cin, line);
        b0_cols = std::stoi(line);
    }
    else {
        current_virus = "1044-1127";
        centralizer_to_check = "D_10";
    }


    std::cout << "Getting T and B0 matrices that worked for group " + centralizer_to_check + " on virus " + current_virus + "...\n";
    std::ifstream fin (curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt");

    if (!fin.is_open()) {
        std::cout << "Failed to open file:\t" + curr_directory + current_virus + "_T_and_B0_pairs_" + centralizer_to_check + ".txt" << std::endl;
        return 0;
    }

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
}