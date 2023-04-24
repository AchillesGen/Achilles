# List of all possible optional binaries in Achilles
option(ACHILLES_ENABLE_CASCADE_TEST "Enable executables for testing the cascade" OFF)
option(ACHILLES_ENABLE_POTENTIAL_TEST "Enable executables for testing the potential" OFF)

# Achilles interface to optional external libraries
option(ACHILLES_ENABLE_SHERPA "Enable the generation of events using Sherpa" OFF) 
option(ACHILLES_ENABLE_ROOT "Enable reading of ROOT flux files" OFF)
SET(ACHILLES_ENABLE_HEPMC3 TRUE)

# BUild options
option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" ON)
option(ACHILLES_ENABLE_GZIP "Enable compression of event files" ON)
option(ACHILLES_ENABLE_PCH "Enable Precompiled Headers" OFF)
option(ACHILLES_ENABLE_TESTING "Enable Test Builds" OFF)
