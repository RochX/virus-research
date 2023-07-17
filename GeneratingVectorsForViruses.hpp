#include <iostream>
#include <fstream>
#include <vector>

#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP
#define VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP

namespace GeneratingVectorsForViruses {
    std::pair<std::vector<int>, std::vector<int>> pickVirusType(const std::string &virus_name, std::vector<EigenType::Vector6f> &starting_generators,
                                                                std::vector<EigenType::Vector6f> &ending_generators, const std::string &curr_directory);

    std::vector<EigenType::Vector6f> generatorsOfPhiX174Native();
    std::vector<EigenType::Vector6f> generatorsOfPhiX174Mature();
    std::vector<EigenType::Vector6f> generatorsOfCCMVNative();
    std::vector<EigenType::Vector6f> generatorsOfCCMVMature();
    std::vector<EigenType::Vector6f> generatorsOfTCVNative();
    std::vector<EigenType::Vector6f> generatorsOfTCVMature();
    std::vector<EigenType::Vector6f> generatorsOfHK97Native();
    std::vector<EigenType::Vector6f> generatorsOfHK97Mature();

    // TODO: implement these
    std::vector<EigenType::Vector6f> generatorsOfCVA10Native();
    std::vector<EigenType::Vector6f> generatorsOfCVA10Mature();
    std::vector<EigenType::Vector6f> generatorsOfCVA10Aparticle();
    std::vector<EigenType::Vector6f> generatorsOfD68Native();
    std::vector<EigenType::Vector6f> generatorsOfD68Mature();
    std::vector<EigenType::Vector6f> generatorsOfD68Aparticle();
    std::vector<EigenType::Vector6f> generatorsOfHE71Native();
    std::vector<EigenType::Vector6f> generatorsOfHE71Mature();
    std::vector<EigenType::Vector6f> generatorsOfHE71Aparticle();

    std::vector<EigenType::Vector6f> startingGeneratorsOfSC_TO_FCC_D10();
    std::vector<EigenType::Vector6f> endingGeneratorsOfSC_TO_FCC_D10();
    std::vector<EigenType::Vector6f> generatorsOf1044();
    std::vector<EigenType::Vector6f> generatorsOf2752();
    std::vector<EigenType::Vector6f> generatorsOf1127();
    std::vector<EigenType::Vector6f> generatorsOf1227();
} // GeneratingVectorsForViruses

#endif //VIRUS_RESEARCH_GENERATINGVECTORSFORVIRUSES_HPP
