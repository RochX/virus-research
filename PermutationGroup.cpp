#include "PermutationGroup.hpp"
#include "Matrix6fFileReader.hpp"

PermutationGroup::PermutationGroup() {
    EigenType::Matrix6f id = EigenType::Matrix6f::Identity();
    EigenType::Matrix6f matrix;

    std::ifstream fin ("permutationGroup.csv");

    if(fin.is_open()) {
        std::cout << "Attempt to read permutation matrices from file..." << std::endl;

        // read in matrices from file
        while(Matrix6fFileReader::readNextMatrix(fin, matrix)) {
            if (std::find(groupElements.begin(), groupElements.end(), matrix) == groupElements.end()) {
                groupElements.push_back(matrix);
            }
        }

        if (groupElements.size() == 720) {
            fin.close();
            return;
        }
    }

    fin.close();


    // this is ugly, but it's the simplest way to generate all permutation matrices...
    std::cout << "No file exists or file reading failed, generating permutation matrices..." << std::endl;
    for (int first = 0; first < 6; ++first) {
        for (int second = 0; second < 6; ++second) {
            for (int third = 0; third < 6; ++third) {
                for (int fourth = 0; fourth < 6; ++fourth) {
                    for (int fifth = 0; fifth < 6; ++fifth) {
                        for (int sixth = 0; sixth < 6; ++sixth) {
                            Eigen::PermutationMatrix<6,6,float> p((Eigen::MatrixXf(6,1) << first, second, third, fourth, fifth, sixth).finished());
                            if (p.toDenseMatrix().colPivHouseholderQr().rank() == 6) {
                                groupElements.push_back(p.toDenseMatrix());
                            }
                        }
                    }
                }
            }
        }
    }

    Eigen::IOFormat CommaSepVals(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
    std::ofstream fout ("permutationGroup.csv");
    for (const EigenType::Matrix6f& P : groupElements) {
        fout << P.format(CommaSepVals) << std::endl;
    }
    fout.close();

    std::cout << "Done." << std::endl;
}

std::string PermutationGroup::groupName() {
    return "PermutationGroup";
}

// NOTE: this does not actually reflect the true centralizer of the permutation group S_6!
bool PermutationGroup::checkIfInCentralizer(EigenType::Matrix6f m) {
    return m.isApprox(EigenType::Matrix6f::Identity());
}