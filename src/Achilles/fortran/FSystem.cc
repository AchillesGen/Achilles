#include "Achilles/System.hh"
#include <iostream>

extern "C" {
char *FindAchillesFile(char *filename, char *head) {
    std::cout << filename << " " << head << std::endl;
    std::string filename_str(filename);
    std::string head_str(head);
    std::string file_loc = achilles::Filesystem::FindFile(filename_str, head_str);

    char *c_file_loc = new char[file_loc.length() + 1];
    std::strcpy(c_file_loc, file_loc.c_str());

    return c_file_loc;
}
}
