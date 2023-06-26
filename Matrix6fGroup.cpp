#include "Matrix6fGroup.hpp"

void Matrix6fGroup::createGroupFromTwoGeneratorsAndKnownSize(const EigenType::Matrix6f &generatorOne, const EigenType::Matrix6f &generatorTwo, int desiredSize) {
    EigenType::Matrix6f currentElement;
    int element_id;
    // generate ICO with a brute force approach as we know it has 60 elements.
    for (int i = 0; i < 100000; ++i) {
        element_id = i;
        currentElement = EigenType::Matrix6f::Identity();

        // generate ICO by using binary representation, with 0 corresponding to A and 1 corresponding to B.
        // exs: 10 = BA
        //      100 = BA^2
        //      1011 = BAB^2
        // note that while we do not get elements like 011 = AB^2 explicitly, we know that B has order 3, and so we do get the element 111011 = B^3AB^2 = AB^2
        // so this process will generate all elements of ICO.
        do {
            if (element_id % 2 == 0) {
                currentElement *= generatorOne;
            }
            else {
                currentElement *= generatorTwo;
            }
            element_id /= 2;
        } while (element_id > 0);

        // check if ICO already has the current element, if not add it.
        if (std::find(groupElements.begin(), groupElements.end(), currentElement) == groupElements.end())
            groupElements.push_back(currentElement);

        // if IcosahedralGroup has 60 (unique) elements we're done, as we know the size of the Icosahedral Group to be 60
        if (groupElements.size() == desiredSize) {
            break;
        }
    }

    if (groupElements.size() != desiredSize) {
        std::cerr << "WARNING: Unable to make group with " << desiredSize << " number of elements. Current size is " << groupElements.size() << std::endl;
    }
}

std::vector<EigenType::Matrix6f> Matrix6fGroup::getGroupElements() {
    return groupElements;
}

std::vector<EigenType::Vector6f> Matrix6fGroup::getOrbitOfVector(const EigenType::Vector6f& vector6f) {
    std::vector<EigenType::Vector6f> orbit;

    for (const EigenType::Matrix6f &matrix : groupElements) {
        // ensure no duplicates in orbit
        if (std::find(orbit.begin(), orbit.end(), matrix*vector6f) == orbit.end())
            orbit.emplace_back(matrix*vector6f);
    }

    return orbit;
}

std::vector<std::vector<EigenType::Vector6f>> Matrix6fGroup::createOrbitsFromGenerators(const std::vector<EigenType::Vector6f>& generators) {
    std::vector<std::vector<EigenType::Vector6f>> orbits;
    for (const EigenType::Vector6f& generator : generators) {
        orbits.push_back(Matrix6fGroup::getOrbitOfVector(generator));
    }
    return orbits;
}