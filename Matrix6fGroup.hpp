#include <vector>
#include "EigenTypes.hpp"

#ifndef VIRUS_RESEARCH_MATRIX6FGROUP_HPP
#define VIRUS_RESEARCH_MATRIX6FGROUP_HPP


class Matrix6fGroup {
protected:
    EigenType::Matrix6f generatorA, generatorB;
    std::vector<EigenType::Matrix6f> groupElements;
    void createGroupFromTwoGeneratorsAndKnownSize(const EigenType::Matrix6f&, const EigenType::Matrix6f&, int);

public:
    std::vector<EigenType::Matrix6f> getGroupElements();
    std::vector<EigenType::Vector6f> getOrbitOfVector(const EigenType::Vector6f&);

    virtual std::string groupName() =0;
    virtual bool checkIfInCentralizer(EigenType::Matrix6f) =0;
};


#endif //VIRUS_RESEARCH_MATRIX6FGROUP_HPP
