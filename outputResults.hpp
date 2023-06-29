#include <iostream>
#include <fstream>

#include "EigenTypes.hpp"
#include "Matrix6fFileReader.hpp"
#include "MatrixFunctions.hpp"

#ifndef VIRUS_RESEARCH_OUTPUTRESULTS_HPP
#define VIRUS_RESEARCH_OUTPUTRESULTS_HPP

namespace outputResults {
    void output3colB0(std::ifstream& fin);
    void output6colB0(std::ifstream& fin);
}

#endif //VIRUS_RESEARCH_OUTPUTRESULTS_HPP
