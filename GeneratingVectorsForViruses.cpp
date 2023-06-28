#include "GeneratingVectorsForViruses.hpp"

namespace GeneratingVectorsForViruses {
    void pickVirusType(const std::string& virus_name, std::vector<EigenType::Vector6f>& starting_generators, std::vector<EigenType::Vector6f>& ending_generators) {
        if (virus_name == "TCV") {
            starting_generators = startingGeneratorsOfTCV();
            ending_generators = endingGeneratorsOfTCV();
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


    std::vector<EigenType::Vector6f> startingGeneratorsOfTCV() {
        EigenType::Matrix6f D;
        D << 1,1,-1,-1,1,1,
                1,1,1,-1,-1,1,
                -1,1,1,1,-1,1,
                -1,-1,1,1,1,1,
                1,-1,-1,1,1,1,
                1,1,1,1,1,1;
        D *= 0.5;

        EigenType::Vector6f s, b, f;
        s << 1, 0, 0, 0, 0, 0;
        b << 1,-1,1,1,-1,1;
        b *= 0.5;
        f << 1,0,0,-1,0,0;
        f *= 0.5;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(b);
        generators.emplace_back(D.inverse() * b);
        generators.emplace_back(D.inverse() * f);

        return generators;
    }

    std::vector<EigenType::Vector6f> endingGeneratorsOfTCV() {
        EigenType::Matrix6f D;
        D << 1,1,-1,-1,1,1,
                1,1,1,-1,-1,1,
                -1,1,1,1,-1,1,
                -1,-1,1,1,1,1,
                1,-1,-1,1,1,1,
                1,1,1,1,1,1;
        D *= 0.5;

        EigenType::Vector6f s, b, f;
        s << 1, 0, 0, 0, 0, 0;
        b << 1,-1,1,1,-1,1;
        b *= 0.5;
        f << 1,0,0,-1,0,0;
        f *= 0.5;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        generators.push_back(b);
        generators.emplace_back(D.inverse() * s);
        generators.emplace_back(D.inverse() * b);
        generators.emplace_back(D.inverse() * f);

        return generators;
    }

    std::vector<EigenType::Vector6f> startingGeneratorsOfSC_TO_FCC_D10() {
        EigenType::Vector6f s;
        s << 1, 0, 0, 0, 0, 0;

        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
        return generators;
    }

    std::vector<EigenType::Vector6f> endingGeneratorsOfSC_TO_FCC_D10() {
        EigenType::Vector6f s, f;
        s << 1, 0, 0, 0, 0, 0;
        f << 1,0,0,-1,0,0;
        f *= 0.5;

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
} // GeneratingVectorsForViruses