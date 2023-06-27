#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_DIHEDRALGROUPONSIX_HPP
#define VIRUS_RESEARCH_DIHEDRALGROUPONSIX_HPP


class DihedralGroupOnSix : public Matrix6fGroup {
public:
    DihedralGroupOnSix();
    bool checkIfInCentralizer(EigenType::Matrix6f) override;
    std::string groupName() override;

    static EigenType::Matrix6f matrixFormOfCentralizer(float x, float y, float z, float t, float u, float w, float v, float s);
};


#endif //VIRUS_RESEARCH_DIHEDRALGROUPONSIX_HPP
