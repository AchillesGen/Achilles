#include "Achilles/System.hh"

extern "C" {

char *FindFile(char *cfilename, char *chead) {
    std::string filename(cfilename);
    std::string head(chead);

    spdlog::info("Looking for file: {} abc", filename);
    std::string result = achilles::Filesystem::FindFile(filename, head);
    spdlog::info("Found file at: {} abc", result);
    char *cresult = new char[result.size()];
    strcpy(cresult, result.c_str());
    return cresult;
}
}
