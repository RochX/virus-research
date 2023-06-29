#include "GeneratingVectorsForViruses.hpp"

EigenType::Matrix6f TwoD {{1,1,-1,-1,1,1},
                       {1, 1, 1, -1, -1, 1},
                       {-1, 1, 1, 1, -1, 1},
                       {-1, -1, 1, 1, 1, 1},
                       {1, -1, -1, 1, 1, 1},
                       {1, 1, 1, 1, 1, 1}};
EigenType::Matrix6f D = 0.5*TwoD;
EigenType::Vector6f s {1, 0, 0, 0, 0, 0};
EigenType::Vector6f b {0.5,-0.5,0.5,0.5, -0.5, 0.5};
EigenType::Vector6f f {0.5,0,0,-0.5,0,0};

namespace GeneratingVectorsForViruses {
    void pickVirusType(const std::string& virus_name, std::vector<EigenType::Vector6f>& starting_generators, std::vector<EigenType::Vector6f>& ending_generators) {
        if (virus_name == "PhiX174") {
            starting_generators = generatorsOfPhiX174Native();
            ending_generators = generatorsOfPhiX174Mature();
        }
        else if (virus_name == "CCMV") {
            starting_generators = generatorsOfCCMVNative();
            ending_generators = generatorsOfCCMVMature();
        }
        else if (virus_name == "TCV") {
            starting_generators = generatorsOfTCVNative();
            ending_generators = generatorsOfTCVMature();
        }
        else if (virus_name == "HK97") {
            starting_generators = generatorsOfHK97Native();
            ending_generators = generatorsOfHK97Mature();
        }
        else if (virus_name == "SC_to_FCC_D10") {
            starting_generators = startingGeneratorsOfSC_TO_FCC_D10();
            ending_generators = endingGeneratorsOfSC_TO_FCC_D10();
        }
        else if (virus_name == "1044-2752") {
            starting_generators = generatorsOf1044();
            ending_generators = generatorsOf2752();
        }
        else if (virus_name == "1044-1127") {
            starting_generators = generatorsOf1044();
            ending_generators = generatorsOf1127();
        }
        else if (virus_name == "1044-1227") {
            starting_generators = generatorsOf1044();
            ending_generators = generatorsOf1227();
        }
    }



    std::vector<EigenType::Vector6f> generatorsOfCCMVNative() {
        std::vector<EigenType::Vector6f> generators;
        generators.emplace_back(D.inverse()*D.inverse()*s);
        generators.emplace_back(b);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfCCMVMature() {
        std::vector<EigenType::Vector6f> generators;
        generators.emplace_back(D*D*b);
        generators.push_back(s);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfTCVNative() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(b);
        generators.emplace_back(D.inverse() * b);
        generators.emplace_back(D.inverse() * f);

        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfTCVMature() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(b);
        generators.emplace_back(D.inverse() * s);
        generators.emplace_back(D.inverse() * b);
        generators.emplace_back(D.inverse() * f);

        return generators;
    }

    std::vector<EigenType::Vector6f> startingGeneratorsOfSC_TO_FCC_D10() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        return generators;
    }

    std::vector<EigenType::Vector6f> endingGeneratorsOfSC_TO_FCC_D10() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(f);
        return generators;
    }
    
    std::vector<EigenType::Vector6f> generatorsOf1044() {
        EigenType::Vector6f v, w, t;
        v << 1, 0, 0, 0, 0, 0;
        w << 1, 1, 0, 0, 1, 1;
        t << 1, 1, 0, 0, 0, 1;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(v);
        generators.push_back(w);
        generators.push_back(t);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOf2752() {
        EigenType::Vector6f v, w, t;
        v << 0, 1, 0, -1, -1, 0;
        w << 0, 0, 1, 0, 1, 0;
        t << 1, 0, 0, 0, 0, 0;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(v);
        generators.push_back(w);
        generators.push_back(t);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOf1127() {
        EigenType::Vector6f v, w, t;
        v << 1, -1, 1, 1, -1, -1;
        v *= 0.5;

        w << 1, 1, -1, 1, 1, -1;
        w *= 0.5;

        t << 3, -1, 1, 1, -1, -1;
        t *= 0.5;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(v);
        generators.push_back(w);
        generators.push_back(t);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOf1227() {
        EigenType::Vector6f v, w, t;

        v << 1, 0, 0, 0, 0, 0;
        w << 1, 1, 0, 0, 0, 1;
        t << 1, 0, 0, 0, 0, 0;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(v);
        generators.push_back(w);
        generators.push_back(t);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfPhiX174Native() {
        std::vector<EigenType::Vector6f> generators;
        generators.emplace_back(D.inverse() * f);
        generators.push_back(s);
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfPhiX174Mature() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(s);
        generators.emplace_back(2*D.inverse() * f);
        // s as translation done
        generators.push_back(b);
        // s as translation done
        generators.push_back(f);
        // s as translation done
        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfHK97Native() {
        std::vector<EigenType::Vector6f> generators;
        generators.emplace_back(D.inverse()*s);
        generators.emplace_back(D.inverse()*f);
        generators.push_back(s);

        return generators;
    }

    std::vector<EigenType::Vector6f> generatorsOfHK97Mature() {
        std::vector<EigenType::Vector6f> generators;
        generators.emplace_back(D.inverse()*s);
        generators.emplace_back(D.inverse()*f);
        generators.emplace_back(D.inverse()*b);
        generators.push_back(s);

        return generators;
    }
} // GeneratingVectorsForViruses