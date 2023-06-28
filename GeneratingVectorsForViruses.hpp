#include <vector>

#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP
#define VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP

namespace GeneratingVectorsForViruses {
    void pickVirusType(const std::string& virus_name, std::vector<EigenType::Vector6f>& starting_generators, std::vector<EigenType::Vector6f>& ending_generators);

    std::vector<EigenType::Vector6f> startingGeneratorsOfTCV();
    std::vector<EigenType::Vector6f> endingGeneratorsOfTCV();
    std::vector<EigenType::Vector6f> startingGeneratorsOfSC_TO_FCC_D10();
    std::vector<EigenType::Vector6f> endingGeneratorsOfSC_TO_FCC_D10();
    std::vector<EigenType::Vector6f> generatorsOf1044();
    std::vector<EigenType::Vector6f> generatorsOf2752();
    std::vector<EigenType::Vector6f> generatorsOf1127();
    std::vector<EigenType::Vector6f> generatorsOf1227();
} // GeneratingVectorsForViruses

#endif //VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP
