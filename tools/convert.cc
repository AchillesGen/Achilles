#include <iostream>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "gzstream/gzstream.h"
#pragma GCC diagnostic pop

#include "H5Cpp.h"
using namespace H5;

std::vector<std::string> tokenize(std::string str, std::string delimiter=" ") {
    size_t pos = 0;
    std::vector<std::string> tokens;
    while((pos = str.find(delimiter)) != std::string::npos) {
        tokens.push_back(str.substr(0, pos));
        str.erase(0, pos + delimiter.size()); 
    }

    return tokens;
}

int main(int argv, char* argc[]) {
    if(argv < 3) {
        std::cout << "Usage: ./convert input output\n";
        std::cout << "   Converts a configuration text file to HDF5\n";
        return -1;
    }
    std::string infile(argc[1]);
    H5std_string outfile(argc[2]);
    const H5std_string dataset_name("Configurations");

    // Load configuration
    igzstream configs(infile.c_str());
    std::string line;
    std::getline(configs, line);
    auto tokens = tokenize(line);
    int nconfigs = std::stoi(tokens[0]);
    double maxWgt = std::stod(tokens[1]);
    double minWgt = std::stod(tokens[2]);

    std::cout << nconfigs << " " << maxWgt << " " << minWgt << std::endl;

    return 0;
}
