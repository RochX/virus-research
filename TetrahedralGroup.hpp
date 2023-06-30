#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP
#define VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP

class TetrahedralGroup: public Matrix6fGroup {
public:
    TetrahedralGroup();
    bool checkIfInCentralizer(EigenType::Matrix6f) override;
    std::string groupName() override;

    static EigenType::Matrix6f matrixFormOfCentralizer(float, float, float, float);
};


#endif //VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP
