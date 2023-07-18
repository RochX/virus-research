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

    std::vector<EigenType::Vector6f> createGeneratorsFromPointArrayNumbers(std::vector<int>& point_array_nums, const std::string& state) {
        std::vector<EigenType::Vector6f> bases;
        std::vector<EigenType::Vector6f> translations;
        std::pair<EigenType::Vector6f, EigenType::Vector6f> base_translation_pair;

        if (!point_array_nums.empty()) {
            std::cout << "The " << state << " state uses config numbers:\t";
            for (int num: point_array_nums)
                std::cout << num << " ";
            std::cout << std::endl;
            int num_to_remove;
            do {
                std::cout << "Enter number to remove, end by entering 0\n";
                std::string line;
                std::getline(std::cin, line);

                try {
                    num_to_remove = std::stoi(line);
                }
                catch (const std::invalid_argument &e) {
                    std::cout << "Invalid number.\n";
                    num_to_remove = -1;
                }
                catch (const std::out_of_range &e) {
                    std::cout << e.what() << std::endl;
                    num_to_remove = -1;
                }

                if (num_to_remove == 0)
                    break;

                auto it = std::find(point_array_nums.begin(), point_array_nums.end(), num_to_remove);
                if (it != point_array_nums.end())
                    point_array_nums.erase(it);
                else
                    std::cout << "Element not found.\n";

            } while (true);
        }
        std::cout << "The " << state << " state will now use numbers:\t";
        for (int num: point_array_nums)
            std::cout << num << " ";
        std::cout << std::endl;

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

    void print_all_configs() {
        initializePointArrayNumbers(point_array_numbering);
        Eigen::IOFormat no_align (Eigen::StreamPrecision, Eigen::DontAlignCols, "\t");
        for (int i = 1; i < point_array_numbering.size(); i++) {
            std::pair<EigenType::Vector6f, EigenType::Vector6f> pair = point_array_numbering[i];
            std::cout << "Config " << i << " base:\t\t" << pair.first.transpose().format(no_align) << std::endl;
//            std::cout << "Config " << i << " translation:\t" << pair.second.transpose().format(no_align) << std::endl;
        }
    }

    void initializeVirusConfigFile(const std::string& VIRUS_CONFIG_FILE_NAME) {
        std::ofstream fout (VIRUS_CONFIG_FILE_NAME);
        fout << "PhiX174\n55\n12,52,29,53\n";
        fout << "CCMV\n10\n27\n";
        fout << "TCV\n8,25,47\n29,13,30,55\n";
        fout << "HK97\n13,55\n13,30,55\n";
        fout << "CVA10_N-M\n1,2\n54,55\n";
        fout << "CVA10_M-A\n54,55\n12,53,13\n";
        fout << "CVA10_N-A\n1,2\n12,53,13\n";
        fout << "D68_N-M\n11,51\n13,54\n";
        fout << "D68_M-A\n13,54\n1,16,3\n";
        fout << "D68_N-A\n11,51\n1,16,3\n";
        fout << "HE71_N-M\n55\n1,3,18\n";
        fout << "HE71_M-A\n1,3,18\n13,54\n";
        fout << "HE71_N-A\n55\n13,54\n";
        fout.close();
    }

    std::vector<int> readIntVector(std::ifstream& fin) {
        std::vector<int> ints;
        std::string line, cell;
        std::getline(fin, line);
        std::stringstream ss (line);

        while (getline(ss, cell, ',')) {
            ints.push_back(std::stoi(cell));
        }

        return ints;
    }
}

namespace GeneratingVectorsForViruses {
    // return a pair indicating what configuration numbers were used
    // returning a value while also modifying generator vectors by reference feels like obfuscation...
    std::pair<std::vector<int>, std::vector<int>> pickVirusType(const std::string &virus_name, std::vector<EigenType::Vector6f> &starting_generators,
                                                                std::vector<EigenType::Vector6f> &ending_generators, const std::string &curr_directory) {
        const std::string VIRUS_CONFIG_FILE_NAME = curr_directory + "virus_configs.txt";
        initializePointArrayNumbers(point_array_numbering);
        initializeVirusConfigFile(VIRUS_CONFIG_FILE_NAME);
        std::ifstream fin (VIRUS_CONFIG_FILE_NAME);
        std::string line;
        std::vector<int> starting_point_array_nums, ending_point_array_nums;
        while (std::getline(fin, line)) {
            if (virus_name == line) {
                starting_point_array_nums = readIntVector(fin);
                ending_point_array_nums = readIntVector(fin);
                break;
            }
        }

        if (!starting_point_array_nums.empty() && !ending_point_array_nums.empty()) {
            starting_generators = createGeneratorsFromPointArrayNumbers(starting_point_array_nums, "starting");
            ending_generators = createGeneratorsFromPointArrayNumbers(ending_point_array_nums, "ending");
            return std::make_pair(starting_point_array_nums, ending_point_array_nums);
        }
        else if (virus_name == "S-F") {
            starting_generators.clear();
            ending_generators.clear();
            starting_generators.push_back(s);
            ending_generators.push_back(f);
        }
        else if (virus_name == "S-B") {
            starting_generators.clear();
            ending_generators.clear();
            starting_generators.push_back(s);
            ending_generators.push_back(b);
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
        else if (virus_name == "print all configs") {
            print_all_configs();
        }

        return std::make_pair(starting_point_array_nums, ending_point_array_nums);
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
} // GeneratingVectorsForViruses