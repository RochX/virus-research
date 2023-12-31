#include "EigenTypes.hpp"
#include "Matrix6fFileReader.hpp"

using std::string;
using std::stringstream;
using std::ifstream;

using namespace EigenType;

// anonymous namespace, these functions are only accessible within this file
namespace {
    // takes a list of lines and converts them to a list of 6x6 matrices
    // precondition: each string in the std::vector lines is 36 comma delimited floats
    std::vector<Matrix6f> convertLinesToMatrices6f(const std::vector<string> &lines) {
        std::vector<Matrix6f> matrices;
        Matrix6f m;

        stringstream ss;
        string cell;
        int pos;
        for (const string &line: lines) {
            ss.clear();
            ss << line;

            pos = 0;
            while (getline(ss, cell, ',')) {
                m((pos / 6) % 6, pos % 6) = stof(cell);
                pos++;
            }

            if (pos != 36) {
                throw pos;
            }

                matrices.push_back(m);
            }


        return matrices;
    }

    std::vector<Matrix6fx3f> convertLinesToMatrices6fx3f(const std::vector<string> &lines) {
        std::vector<Matrix6fx3f> matrices;
        Matrix6fx3f m;

        stringstream ss;
        string cell;
        int pos;
        for (const string &line: lines) {
            ss.clear();
            ss << line;

            pos = 0;
            while (getline(ss, cell, ',')) {
                m((pos / 3) % 6, pos % 3) = stof(cell);
                pos++;
            }

            if (pos != 18) {
                throw pos;
            }

            matrices.push_back(m);
        }


        return matrices;
    }

    // reads the next N lines from the specified ifstream
    // the ifstream has less than N lines left, it will read the rest of them
    bool readNextNLines(ifstream &fin, std::vector<string> &lines, int N) {
        string line;
        int lines_read = 0;
        while (getline(fin, line)) {
            lines.push_back(line);

            lines_read++;
            if (lines_read >= N) {
                break;
            }
        }

        return lines_read > 0;
    }
} // end anonymous namespace

namespace Matrix6fFileReader {
    bool readNextMatrix(ifstream &fin, Matrix6f &m) {
        std::vector<Matrix6f> matrices;

        // try to read matrix
        if (readNextNMatrices(fin, matrices, 1)) {
            // succeeded, assign it
            m = matrices.front();
            return true;
        }

        // failed to read matrix
        return false;
    }

    bool readNextNMatrices(ifstream &fin, std::vector<Matrix6f> &matrices, int N) {
        std::vector<string> lines;

        readNextNLines(fin, lines, N);

        try {
            matrices = convertLinesToMatrices6f(lines);
        }
        catch(int numElements) {
            std::cerr << "Incorrect number of elements in Matrix6fReader.\nGot " << numElements << " when expecting 36." << std::endl;
            return false;
        }

        return !matrices.empty();
    }
}

namespace Matrix6fx3fFileReader {
    bool readNextMatrix(std::ifstream &fin, EigenType::Matrix6fx3f &m) {
        std::vector<Matrix6fx3f> matrices;

        // try to read matrix
        if (readNextNMatrices(fin, matrices, 1)) {
            // succeeded, assign it
            m = matrices.front();
            return true;
        }

        // failed to read matrix
        return false;
    }

    bool readNextNMatrices(std::ifstream &fin, std::vector<EigenType::Matrix6fx3f> &matrices, int N) {
        std::vector<string> lines;

        readNextNLines(fin, lines, N);

        try {
            matrices = convertLinesToMatrices6fx3f(lines);
        }
        catch(int numElements) {
            std::cerr << "Incorrect number of elements in Matrix6fx3fReader.\nGot " << numElements << " when expecting 18." << std::endl;
            return false;
        }

        return !matrices.empty();
    }
}



