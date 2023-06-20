#include <iostream>
#include <fstream>
#include <vector>
#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_MATRIX6FFILEREADER_H
#define VIRUS_RESEARCH_MATRIX6FFILEREADER_H

namespace Matrix6fFileReader {
    bool readNextMatrix(std::ifstream &, EigenType::Matrix6f &);
    bool readNextNMatrices(std::ifstream &, std::vector<EigenType::Matrix6f> &, int);
}

#endif //VIRUS_RESEARCH_MATRIX6FFILEREADER_H
