#include "GeneratingVectorsForViruses.hpp"

std::vector<std::pair<EigenType::Vector6f, EigenType::Vector6f>> point_array_numbering;
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

namespace {
    void initializePointArrayNumbers(std::vector<std::pair<EigenType::Vector6f, EigenType::Vector6f>>& point_array_numbering) {
        point_array_numbering.resize(56);
        point_array_numbering[0] = std::make_pair(EigenType::Vector6f::Zero(), EigenType::Vector6f::Zero()); // use as useless assignment to not be off by one
        point_array_numbering[1] = std::make_pair(D*s, f);
        point_array_numbering[2] = std::make_pair(0.5*D*D*s, f);
        point_array_numbering[3] = std::make_pair(s, f);
        point_array_numbering[4] = std::make_pair(0.5*D*s, f);
        point_array_numbering[5] = std::make_pair(0.5*s, f);
        point_array_numbering[6] = std::make_pair(0.5*D.inverse()*s, f);
        point_array_numbering[7] = std::make_pair(D*s, b);
        point_array_numbering[8] = std::make_pair(s, b);
        point_array_numbering[9] = std::make_pair(D.inverse()*s, b);
        point_array_numbering[10] = std::make_pair(D.inverse()*D.inverse()*s, b);
        point_array_numbering[11] = std::make_pair(D*s, s);
        point_array_numbering[12] = std::make_pair(s, s);
        point_array_numbering[13] = std::make_pair(D.inverse()*s, s);
        point_array_numbering[14] = std::make_pair(D*D*b, f);
        point_array_numbering[15] = std::make_pair(0.5*D*D*D*b, f);
        point_array_numbering[16] = std::make_pair(D*b, f);
        point_array_numbering[17] = std::make_pair(0.5*D*D*b, f);
        point_array_numbering[18] = std::make_pair(b, f);
        point_array_numbering[19] = std::make_pair(0.5*D*b, f);
        point_array_numbering[20] = std::make_pair(0.5*b, f);
        point_array_numbering[21] = std::make_pair(0.5*D.inverse()*b, f);
        point_array_numbering[22] = std::make_pair(D*D*b, b);
        point_array_numbering[23] = std::make_pair(D*b, b);
        point_array_numbering[24] = std::make_pair(b, b);
        point_array_numbering[25] = std::make_pair(D.inverse()*b, b);
        point_array_numbering[26] = std::make_pair(D.inverse()*D.inverse()*b, b);
        point_array_numbering[27] = std::make_pair(D*D*b, s);
        point_array_numbering[28] = std::make_pair(D*b, s);
        point_array_numbering[29] = std::make_pair(b, s);
        point_array_numbering[30] = std::make_pair(D.inverse()*b, s);
        point_array_numbering[31] = std::make_pair(2*D*f, f);
        point_array_numbering[32] = std::make_pair(D*D*f, f);
        point_array_numbering[33] = std::make_pair(2*f, f);
        point_array_numbering[34] = std::make_pair(D*f, f);
        point_array_numbering[35] = std::make_pair(2*D.inverse()*f, f);
        point_array_numbering[36] = std::make_pair(f, f);
        point_array_numbering[37] = std::make_pair(0.5*D*f, f);
        point_array_numbering[38] = std::make_pair(D.inverse()*f, f);
        point_array_numbering[39] = std::make_pair(0.5*f, f);
        point_array_numbering[40] = std::make_pair(D.inverse()*D.inverse()*f, f);
        point_array_numbering[41] = std::make_pair(0.5*D.inverse()*f, f);
        point_array_numbering[42] = std::make_pair(2*D*f, b);
        point_array_numbering[43] = std::make_pair(2*f, b);
        point_array_numbering[44] = std::make_pair(2*D.inverse()*f, b);
        point_array_numbering[45] = std::make_pair(f, b);
        point_array_numbering[46] = std::make_pair(2*D.inverse()*D.inverse()*f, b);
        point_array_numbering[47] = std::make_pair(D.inverse()*f, b);
        point_array_numbering[48] = std::make_pair(2*D.inverse()*D.inverse()*D.inverse()*f, b);
        point_array_numbering[49] = std::make_pair(D.inverse()*D.inverse()*f, b);
        point_array_numbering[50] = std::make_pair(2*D*f, s);
        point_array_numbering[51] = std::make_pair(2*f, s);
        point_array_numbering[52] = std::make_pair(2*D.inverse()*f, s);
        point_array_numbering[53] = std::make_pair(f, s);
        point_array_numbering[54] = std::make_pair(2*D.inverse()*D.inverse()*f, s);
        point_array_numbering[55] = std::make_pair(D.inverse()*f, s);
    }

    std::vector<EigenType::Vector6f> createGeneratorsFromPointArrayNumbers(const std::vector<int>& point_array_nums) {
        std::vector<EigenType::Vector6f> bases;
        std::vector<EigenType::Vector6f> translations;
        std::pair<EigenType::Vector6f, EigenType::Vector6f> base_translation_pair;
        for (int point_array_num : point_array_nums) {
            base_translation_pair = point_array_numbering[point_array_num];

            // we can have multiple of the same base
            bases.push_back(base_translation_pair.first);

            // we cannot have multiple of the same translation
            if (std::find(translations.begin(), translations.end(), base_translation_pair.second) == translations.end()) {
                translations.push_back(base_translation_pair.second);
            }
        }
        // should only be one translation
        assert(translations.size() == 1);

        // convention: put translation vector first!
        std::vector<EigenType::Vector6f> generators;
        generators.insert(generators.end(), translations.begin(), translations.end());
        generators.insert(generators.end(), bases.begin(), bases.end());

        return generators;
    }
}

namespace GeneratingVectorsForViruses {
    void pickVirusType(const std::string& virus_name, std::vector<EigenType::Vector6f>& starting_generators, std::vector<EigenType::Vector6f>& ending_generators) {
        initializePointArrayNumbers(point_array_numbering);
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
        else if (virus_name == "SC_FCC_D10") {
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
        else if (virus_name == "CVA10_N-M") {
            starting_generators = generatorsOfCVA10Native();
            ending_generators = generatorsOfCVA10Mature();
        }
        else if (virus_name == "CVA10_M-A") {
            starting_generators = generatorsOfCVA10Mature();
            ending_generators = generatorsOfCVA10Aparticle();
        }
        else if (virus_name == "CVA10_N-A") {
            starting_generators = generatorsOfCVA10Native();
            ending_generators = generatorsOfCVA10Aparticle();
        }
        else if (virus_name == "D68_N-M") {
            starting_generators = generatorsOfD68Native();
            ending_generators = generatorsOfD68Mature();
        }
        else if (virus_name == "D68_M-A") {
            starting_generators = generatorsOfD68Mature();
            ending_generators = generatorsOfD68Aparticle();
        }
        else if (virus_name == "D68_N-A") {
            starting_generators = generatorsOfD68Native();
            ending_generators = generatorsOfD68Aparticle();
        }
        else if (virus_name == "HE71_N-M") {
            starting_generators = generatorsOfHE71Native();
            ending_generators = generatorsOfHE71Mature();
        }
        else if (virus_name == "HE71_M-A") {
            starting_generators = generatorsOfHE71Mature();
            ending_generators = generatorsOfHE71Aparticle();
        }
        else if (virus_name == "HE71_N-A") {
            starting_generators = generatorsOfHE71Native();
            ending_generators = generatorsOfHE71Aparticle();
        }
    }



    std::vector<EigenType::Vector6f> generatorsOfCCMVNative() {
        std::vector<int> point_array_nums {10};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfCCMVMature() {
        std::vector<int> point_array_nums {27};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfTCVNative() {
        std::vector<int> point_array_nums {8, 25, 47};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfTCVMature() {
        std::vector<int> point_array_nums {29, 13, 30, 55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> startingGeneratorsOfSC_TO_FCC_D10() {
        std::vector<EigenType::Vector6f> generators;
        generators.push_back(s);
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
        std::vector<int> point_array_nums {55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfPhiX174Mature() {
        std::vector<int> point_array_nums {12, 52, 29, 53};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfHK97Native() {
        std::vector<int> point_array_nums {13, 55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfHK97Mature() {
        std::vector<int> point_array_nums {13, 30, 55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfCVA10Native() {
        std::vector<int> point_array_nums {1, 2};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
    std::vector<EigenType::Vector6f> generatorsOfCVA10Mature() {
        std::vector<int> point_array_nums {54, 55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
    std::vector<EigenType::Vector6f> generatorsOfCVA10Aparticle() {
        std::vector<int> point_array_nums {12, 53, 13};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfD68Native() {
        std::vector<int> point_array_nums {11, 51};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
    std::vector<EigenType::Vector6f> generatorsOfD68Mature() {
        std::vector<int> point_array_nums {13, 54};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
    std::vector<EigenType::Vector6f> generatorsOfD68Aparticle() {
        std::vector<int> point_array_nums {1, 16, 3};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfHE71Native() {
        std::vector<int> point_array_nums {55};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }

    std::vector<EigenType::Vector6f> generatorsOfHE71Mature() {
        std::vector<int> point_array_nums {1, 3, 18};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
    std::vector<EigenType::Vector6f> generatorsOfHE71Aparticle() {
        std::vector<int> point_array_nums {13, 54};
        return createGeneratorsFromPointArrayNumbers(point_array_nums);
    }
} // GeneratingVectorsForViruses