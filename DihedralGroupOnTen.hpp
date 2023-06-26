#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_DIHEDRALGROUPONTEN_HPP
#define VIRUS_RESEARCH_DIHEDRALGROUPONTEN_HPP


class DihedralGroupOnTen : public Matrix6fGroup {
public:
    DihedralGroupOnTen();
    bool checkIfInCentralizer(EigenType::Matrix6f) override;
    std::string groupName() override;

    static EigenType::Matrix6f matrixFormOfCentralizer(float x, float y, float z, float t, float u, float w);
};


#endif //VIRUS_RESEARCH_DIHEDRALGROUPONTEN_HPP
