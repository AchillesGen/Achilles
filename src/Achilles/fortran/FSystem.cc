#include "Achilles/System.hh"

extern "C" {

char *FindFile(char *cfilename, char *chead) {
    std::string filename(cfilename);
    std::string head(chead);

    std::string result = achilles::Filesystem::FindFile(filename, head);
    char *cresult = new char[result.size()];
    strcpy(cresult, result.c_str());
    return cresult;
}
}
