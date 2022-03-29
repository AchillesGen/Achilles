#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

#include <iostream>

using namespace HepMC3;

int main(int argc, char **argv) {
    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " <HepMC3_input_file>" << std::endl;
        return -1;
    }

    ReaderAscii input_file(argv[1]);

    size_t parsed_events = 0;

    while(!input_file.failed()) {
        GenEvent evt(Units::MEV, Units::MM);

        // Read event from input file
        input_file.read_event(evt);
        if(input_file.failed()) break;

        std::cout << "Event #: " << ++parsed_events <<  " Event wgt: " << evt.weight() << std::endl;
        for(auto p : evt.particles()) {
            auto mom = p -> momentum();
            std::cout << " " << p -> pid() << " " << p -> status() << " " << mom.px() << " " << mom.py() << " " << mom.pz() << " " << mom.e() << std::endl;
        }
    }

    input_file.close();
    return 0;
}
