#include "Achilles/ReferenceHandler.hh"
#include "fmt/format.h"
#include <fstream>

using achilles::ReferenceHandler;

void ReferenceHandler::WriteReferences() const {
    fmt::print("Writing out references from current run.\n");
    fmt::print("References can be found in \"References.txt\".\n");
    std::ofstream outfile("References.txt");
    outfile << "Reference List:\n";
    for(const auto &ref : m_references) {
        outfile << fmt::format("- ID: {}\n  Description: {}\n", ref.second.id,
                               ref.second.description);
    }
}
