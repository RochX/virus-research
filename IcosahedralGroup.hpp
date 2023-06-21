#include <vector>
#include "EigenTypes.hpp"
#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
#define VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP

class IcosahedralGroup : public Matrix6fGroup {
public:
    IcosahedralGroup();
    bool checkIfInCentralizer(EigenType::Matrix6f) override;
    std::string groupName() override;

    static EigenType::Matrix6f matrixFormOfCentralizer(float, float);
};


#endif //VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
