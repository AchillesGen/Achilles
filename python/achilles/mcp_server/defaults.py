# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later

"""Constants and mappings for the Achilles MCP server."""

# Beam type -> PDG ID and process lepton configuration
BEAM_TYPES = {
    "electron": {"pid": 11, "leptons": [11, [11]]},
    "nue": {"pid": 12, "leptons": [12, [11]]},
    "nuebar": {"pid": -12, "leptons": [-12, [-11]]},
    "numu": {"pid": 14, "leptons": [14, [13]]},
    "numubar": {"pid": -14, "leptons": [-14, [-13]]},
}

AVAILABLE_NUCLEI = ["12C", "16O", "40Ar", "1H", "1N"]

NUCLEAR_MODELS = ["QESpectral", "FortranModel"]

FORTRAN_MODEL_NAMES = [
    "QE_Spectral_Func",
    "RES_Spectral_Func",
    "Intf_Spectral_Func",
]

FORM_FACTOR_VECTOR = ["Kelly", "BBBA", "ArringtonHill", "VectorDipole"]
FORM_FACTOR_AXIAL = ["AxialDipole", "AxialZExpansion"]
FORM_FACTOR_COHERENT = ["Helm", "Lovato"]

OUTPUT_FORMATS = ["NuHepMC", "HepMC3"]

# Beam parameter types
# - Monochromatic: fixed energy, requires "Energy" (MeV)
# - Spectrum: energy distribution from file, requires one of:
#     "Histogram" (path to .dat file) or "HepData" (path to .yaml file)
# - FlatFlux: uniform distribution, requires "MinEnergy" and "MaxEnergy" (MeV)
BEAM_PARAM_TYPES = ["Monochromatic", "Spectrum", "FlatFlux"]

CASCADE_PROBABILITY_TYPES = ["Cylinder", "Gaussian"]

# Nuclei that have spectral function data files available
NUCLEUS_SPECTRAL_FUNCTIONS = {
    "12C": {
        "SpectralP": "data/Spectral_Functions/pke12p_tot.data",
        "SpectralN": "data/Spectral_Functions/pke12n_tot.data",
    },
    "16O": {
        "SpectralP": "data/Spectral_Functions/pke16p_tot.data",
        "SpectralN": "data/Spectral_Functions/pke16n_tot.data",
    },
    "40Ar": {
        "SpectralP": "data/Spectral_Functions/pke40p_tot.data",
        "SpectralN": "data/Spectral_Functions/pke40n_tot.data",
    },
    "1H": None,
    "1N": None,
}

SINGLE_NUCLEON_TARGETS = {"1H", "1N"}

VALID_CUT_TYPES = {"AngleTheta", "Energy", "Momentum", "TransverseMomentum", "ETheta2"}
