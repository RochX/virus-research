#include <vector>
#include "EigenTypes.hpp"
#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
#define VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP

class IcosahedralGroup : public Matrix6fGroup {
private:
    static EigenType::Matrix6f matrixFormOfCentralizer(float, float);

public:
    IcosahedralGroup();
    static bool checkIfInCentralizer(EigenType::Matrix6f);
};


#endif //VIRUS_RESEARCH_ICOSAHEDRALGROUP_HPP
