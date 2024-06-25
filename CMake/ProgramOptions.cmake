# List of all possible optional binaries in Achilles
option(ACHILLES_ENABLE_CASCADE_TEST "Enable executables for testing the cascade" OFF)
option(ACHILLES_ENABLE_POTENTIAL_TEST "Enable executables for testing the potential" OFF)
option(ACHILLES_ENABLE_PRECOMPUTED "Enable executables for running the cascade on precomputed events" OFF)

# Achilles interface to optional external libraries
option(ACHILLES_ENABLE_SHERPA "Enable the generation of events using Sherpa" OFF) 
option(ACHILLES_ENABLE_ROOT "Enable reading of ROOT flux files" OFF)
SET(ACHILLES_ENABLE_HEPMC3 TRUE)

# BUild options
option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)
option(ACHILLES_ENABLE_GZIP "Enable compression of event files" ON)
option(ACHILLES_ENABLE_PCH "Enable Precompiled Headers" OFF)
option(ACHILLES_ENABLE_TESTING "Enable Test Builds" OFF)

# Achilles options for debugging
option(ACHILLES_LOW_MEMORY "Reduce Achilles memory usage at cost of performance" OFF)
option(ACHILLES_EVENT_DETAILS "Produces all event details when in trace mode" OFF)
