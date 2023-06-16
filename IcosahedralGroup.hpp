#include <vector>
#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
#define VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP

class IcosahedralGroup {
private:
    EigenType::Matrix6f generatorA, generatorB;
    std::vector<EigenType::Matrix6f> groupElements;
    static EigenType::Matrix6f matrixFormOfCentralizer(float, float);

public:
    IcosahedralGroup();
    std::vector<EigenType::Matrix6f> getGroupElements();
    std::vector<EigenType::Vector6f> getOrbitOfVector(EigenType::Vector6f);
    static bool checkIfInCentralizer(EigenType::Matrix6f);
};


#endif //VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
