#include "GeneratingVectorsForViruses.hpp"

namespace GeneratingVectorsForViruses {
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
} // GeneratingVectorsForViruses