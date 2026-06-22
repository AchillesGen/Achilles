"""Run card builder with YAML serialization supporting !include directives."""

import yaml

from .defaults import BEAM_TYPES, AVAILABLE_NUCLEI


class IncludeTag:
    """Sentinel for !include directives in YAML output."""

    def __init__(self, path: str):
        self.path = path

    def __repr__(self):
        return f"!include \"{self.path}\""


class NoneString:
    """Sentinel to emit the string 'None' instead of YAML null."""

    def __repr__(self):
        return "None"


def _include_representer(dumper, data):
    return dumper.represent_scalar("!include", data.path)


def _none_string_representer(dumper, data):
    return dumper.represent_scalar("tag:yaml.org,2002:str", "None")


class AchillesDumper(yaml.SafeDumper):
    pass


AchillesDumper.add_representer(IncludeTag, _include_representer)
AchillesDumper.add_representer(NoneString, _none_string_representer)


class RunCardBuilder:
    """Builds Achilles run card configurations as Python dicts, serializable to YAML."""

    def __init__(self):
        self._config: dict = {}
        self.reset()

    def reset(self):
        """Reset to sensible defaults (nue on 12C, 1000 MeV monochromatic)."""
        self._config = {
            "Main": {
                "NEvents": 10000,
                "HardCuts": False,
                "EventCuts": False,
                "DoRotate": False,
                "RunDecays": False,
                "Output": {
                    "Format": "NuHepMC",
                    "Name": "achilles.hepmc",
                    "Zipped": True,
                },
            },
            "SherpaOptions": IncludeTag("data/default/SherpaOptions.yml"),
            "Processes": [{"Leptons": [12, [11]]}],
            "Beams": [
                {
                    "Beam": {
                        "PID": 12,
                        "Beam Params": {
                            "Type": "Monochromatic",
                            "Energy": 1000,
                        },
                    }
                }
            ],
            "Cascade": {
                "Run": False,
                "Interactions": IncludeTag(
                    "data/default/VirtResInteractions.yml"
                ),
                "Step": 0.04,
                "Probability": "Cylinder",
                "InMedium": NoneString(),
                "PotentialProp": False,
                "Algorithm": "Base",
            },
            "NuclearModels": [
                {
                    "NuclearModel": {
                        "Model": "QESpectral",
                        "FormFactorFile": "FormFactors.yml",
                        "SpectralP": "data/Spectral_Functions/pke12p_tot.data",
                        "SpectralN": "data/Spectral_Functions/pke12n_tot.data",
                        "Ward": NoneString(),
                    }
                }
            ],
            "Nuclei": [{"Nucleus": IncludeTag("data/default/12C.yml")}],
            "HardCuts": [],
            "Options": {
                "Initialize": {"Seed": 12345678, "Accuracy": 1e-2},
                "Unweighting": {"Name": "Percentile", "percentile": 99},
            },
            "Backend": {"Name": "Default", "Options": []},
        }

    def set_main(
        self,
        nevents: int | None = None,
        hard_cuts_enabled: bool | None = None,
        event_cuts_enabled: bool | None = None,
        do_rotate: bool | None = None,
        run_decays: bool | None = None,
        output_format: str | None = None,
        output_name: str | None = None,
        output_zipped: bool | None = None,
    ):
        main = self._config["Main"]
        if nevents is not None:
            main["NEvents"] = nevents
        if hard_cuts_enabled is not None:
            main["HardCuts"] = hard_cuts_enabled
        if event_cuts_enabled is not None:
            main["EventCuts"] = event_cuts_enabled
        if do_rotate is not None:
            main["DoRotate"] = do_rotate
        if run_decays is not None:
            main["RunDecays"] = run_decays
        output = main["Output"]
        if output_format is not None:
            output["Format"] = output_format
        if output_name is not None:
            output["Name"] = output_name
        if output_zipped is not None:
            output["Zipped"] = output_zipped

    def set_beam(
        self,
        beam_type: str,
        energy: float | None = None,
        beam_params_type: str = "Monochromatic",
        histogram: str | None = None,
        hepdata: str | None = None,
        min_energy: float | None = None,
        max_energy: float | None = None,
    ):
        """Set beam configuration.

        For Monochromatic: provide energy.
        For Spectrum: provide histogram (path to .dat file) or hepdata (path to .yaml file).
        For FlatFlux: provide min_energy and max_energy.
        """
        if beam_type not in BEAM_TYPES:
            raise ValueError(
                f"Unknown beam_type '{beam_type}'. "
                f"Choose from: {list(BEAM_TYPES.keys())}"
            )
        info = BEAM_TYPES[beam_type]

        beam_params: dict = {"Type": beam_params_type}
        if beam_params_type == "Monochromatic":
            if energy is None:
                raise ValueError("energy is required for Monochromatic beams")
            beam_params["Energy"] = energy
        elif beam_params_type == "Spectrum":
            if histogram:
                beam_params["Histogram"] = histogram
            elif hepdata:
                beam_params["HepData"] = hepdata
            else:
                raise ValueError(
                    "Spectrum beams require either histogram or hepdata path"
                )
        elif beam_params_type == "FlatFlux":
            if min_energy is None or max_energy is None:
                raise ValueError(
                    "FlatFlux beams require both min_energy and max_energy"
                )
            beam_params["MinEnergy"] = min_energy
            beam_params["MaxEnergy"] = max_energy

        self._config["Beams"] = [
            {
                "Beam": {
                    "PID": info["pid"],
                    "Beam Params": beam_params,
                }
            }
        ]
        self._config["Processes"] = [{"Leptons": info["leptons"]}]

    def set_nucleus(self, name: str, use_defaults: bool = True):
        if name not in AVAILABLE_NUCLEI:
            raise ValueError(
                f"Unknown nucleus '{name}'. "
                f"Choose from: {AVAILABLE_NUCLEI}"
            )
        if use_defaults:
            self._config["Nuclei"] = [
                {"Nucleus": IncludeTag(f"data/default/{name}.yml")}
            ]
        else:
            self._config["Nuclei"] = [{"Nucleus": {"Name": name}}]

    def set_nuclear_model(
        self,
        model: str,
        fortran_model_name: str | None = None,
        form_factor_file: str = "FormFactors.yml",
        spectral_p: str | None = None,
        spectral_n: str | None = None,
        config_file: str | None = None,
        ward: str = "None",
    ):
        """Set the nuclear model, replacing any existing models."""
        model_config = self._build_nuclear_model(
            model, fortran_model_name, form_factor_file,
            spectral_p, spectral_n, config_file, ward,
        )
        self._config["NuclearModels"] = [{"NuclearModel": model_config}]

    def add_nuclear_model(
        self,
        model: str,
        fortran_model_name: str | None = None,
        form_factor_file: str = "FormFactors.yml",
        spectral_p: str | None = None,
        spectral_n: str | None = None,
        config_file: str | None = None,
        ward: str = "None",
    ):
        """Append an additional nuclear model."""
        model_config = self._build_nuclear_model(
            model, fortran_model_name, form_factor_file,
            spectral_p, spectral_n, config_file, ward,
        )
        self._config["NuclearModels"].append({"NuclearModel": model_config})

    def _build_nuclear_model(
        self, model, fortran_model_name, form_factor_file,
        spectral_p, spectral_n, config_file, ward,
    ):
        config = {"Model": model, "FormFactorFile": form_factor_file}
        if model == "FortranModel":
            if not fortran_model_name:
                raise ValueError(
                    "fortran_model_name is required when model='FortranModel'"
                )
            config["Name"] = fortran_model_name
        if spectral_p:
            config["SpectralP"] = spectral_p
        if spectral_n:
            config["SpectralN"] = spectral_n
        if config_file:
            config["ConfigFile"] = config_file
        config["Ward"] = NoneString() if ward == "None" else ward
        return config

    def set_cascade(
        self,
        run: bool,
        step: float = 0.04,
        probability: str = "Cylinder",
        in_medium: str = "None",
        potential_prop: bool = False,
        algorithm: str = "Base",
    ):
        self._config["Cascade"] = {
            "Run": run,
            "Interactions": IncludeTag(
                "data/default/VirtResInteractions.yml"
            ),
            "Step": step,
            "Probability": probability,
            "InMedium": NoneString() if in_medium == "None" else in_medium,
            "PotentialProp": potential_prop,
            "Algorithm": algorithm,
        }

    def set_hard_cuts(self, cuts: list[dict]):
        self._config["HardCuts"] = cuts
        self._config["Main"]["HardCuts"] = len(cuts) > 0

    def set_initialize(
        self,
        seed: int | None = None,
        accuracy: float | None = None,
    ):
        init = self._config["Options"]["Initialize"]
        if seed is not None:
            init["Seed"] = seed
        if accuracy is not None:
            init["Accuracy"] = accuracy

    def set_unweighting(
        self,
        name: str = "Percentile",
        percentile: float = 99,
    ):
        self._config["Options"]["Unweighting"] = {
            "Name": name,
            "percentile": percentile,
        }

    def get_config(self) -> dict:
        return self._config

    def to_yaml(self) -> str:
        return yaml.dump(
            self._config,
            Dumper=AchillesDumper,
            default_flow_style=False,
            sort_keys=False,
        )
