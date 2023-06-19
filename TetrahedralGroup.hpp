#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP
#define VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP

class TetrahedralGroup: public Matrix6fGroup {
private:
    static EigenType::Matrix6f matrixFormOfCentralizer(float, float, float, float);

public:
    TetrahedralGroup();
    static bool checkIfInCentralizer(EigenType::Matrix6f);
};


#endif //VIRUS_RESEARCH_TETRAHEDRALGROUP_HPP
