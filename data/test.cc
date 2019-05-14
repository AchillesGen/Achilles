#include <iostream>
#include <vector>

#include "GeantPPData.hh"
#include "GeantNPData.hh"

int main() {
    for(size_t i = 0; i < G4LEpp::Sig.size(); ++i) {
        for(size_t j = 0; j < G4LEpp::Sig[i].size(); ++j) {
            std::cout << G4LEpp::Sig[i][j]/G4LEpp::Sigtot[i] << ", ";
        }
        std::cout << std::endl;
    }

    return 0;
}
