#include "Achilles/Debug.hh"
#include "Achilles/Utilities.hh"
#include <fstream>

achilles::DebugEvents::DebugEvents(const std::string &filename, size_t multiplicity) {
    std::ifstream input(filename.c_str());
    std::vector<std::string> tokens;
    std::string line;
    size_t nline = 0, index;
    std::vector<FourVector> mom(multiplicity);
    while(std::getline(input, line)) {
        index = 0;
        nline++;
        tokens.clear();
        tokenize(line, tokens, ", ");
        if(tokens.size() % 4 != 0 || tokens.size() != 4 * multiplicity) {
            auto msg = fmt::format("DebugEvents: Invalid event at line {}", nline);
            throw std::runtime_error(msg);
        }
        for(auto &p : mom) {
            for(size_t i = 0; i < 4; ++i) { p[i] = std::stod(tokens[index++]); }
        }
        events.push_back(mom);
    }
}
