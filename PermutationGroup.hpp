#include <fstream>
#include "Matrix6fGroup.hpp"

#ifndef VIRUS_RESEARCH_PERMUTATIONGROUP_HPP
#define VIRUS_RESEARCH_PERMUTATIONGROUP_HPP

class PermutationGroup : public Matrix6fGroup {
public:
    PermutationGroup();
    std::string groupName() override;

    bool checkIfInCentralizer(EigenType::Matrix6f m);
};

#endif //VIRUS_RESEARCH_PERMUTATIONGROUP_HPP
