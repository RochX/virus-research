#include "outputResults.hpp"

using namespace EigenType;

void outputResults::output3colB0(std::ifstream& fin) {
    Matrix6f transition;
    Matrix6fx3f b0, prod;
    std::string count_s, blank;

    while (getline(fin, count_s)) {
        Matrix6fFileReader::readNextMatrix(fin, transition);
        Matrix6fx3fFileReader::readNextMatrix(fin, b0);
        Matrix6fx3fFileReader::readNextMatrix(fin, prod);
        getline(fin, blank); // read blank line

        MatrixFunctions::fixZeroEntries(b0);
        MatrixFunctions::fixZeroEntries(prod);

        std::cout << count_s << ":" << std::endl;
        std::cout << "\tT:" << std::endl;
        std::cout << transition.format(EigenType::TAB_INDENT) << std::endl;
        std::cout << "\tB0:" << std::endl;
        std::cout << b0.format(EigenType::TAB_INDENT) << std::endl;
        std::cout << "\tT*B0:" << std::endl;
        std::cout << prod.format(EigenType::TAB_INDENT) << std::endl << std::endl;
    }
    std::cout << "Done with file.\n";
    fin.close();
}

void outputResults::output6colB0(std::ifstream& fin) {
    Matrix6f transition;
    Matrix6f b0, prod;
    std::string count_s, blank;

    while (getline(fin, count_s)) {
        Matrix6fFileReader::readNextMatrix(fin, transition);
        Matrix6fFileReader::readNextMatrix(fin, b0);
        Matrix6fFileReader::readNextMatrix(fin, prod);
        getline(fin, blank); // read blank line

        MatrixFunctions::fixZeroEntries(b0);
        MatrixFunctions::fixZeroEntries(prod);

        std::cout << count_s << ":" << std::endl;
        std::cout << "\tT:" << std::endl;
        std::cout << transition.format(EigenType::TAB_INDENT) << std::endl;
        std::cout << "\tB0:" << std::endl;
        std::cout << b0.format(EigenType::TAB_INDENT) << std::endl;
        std::cout << "\tT*B0:" << std::endl;
        std::cout << prod.format(EigenType::TAB_INDENT) << std::endl << std::endl;
    }
    std::cout << "Done with file.\n";
    fin.close();
}