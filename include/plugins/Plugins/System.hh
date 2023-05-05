#ifndef SYSTEMVARIABLES_HH
#define SYSTEMVARIABLES_HH

#include <string>
#if __cplusplus == 201402L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else // if __cplusplus == 201703L
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace achilles {

namespace SystemVariables {
const std::string libraryPrefix = "lib";
const std::string librarySuffix = ".so";
} // namespace SystemVariables

class Directory {};

} // namespace achilles

#endif
