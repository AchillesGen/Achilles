# Sudakov veto cascade: validation

Code: `include/Achilles/Sudakov.hh`, `src/Achilles/Sudakov.cc` (engine, Track A `DensityRate`,
Track B `ConfigurationRate`, `ForwardCollision`, `AchillesCollision`), Cascade algorithms
`SudakovA` / `SudakovB` (`Cascade::EvolveSudakov`), unit tests `test/test_sudakov.cc`.

Build with `-DACHILLES_ENABLE_CASCADE_TEST=ON` to get `achilles-sudakov-validate`.

1. Step 4(a)-(e), forward (no-deflection) mode against closed forms, exact enumeration and the
   QMC reference numbers:
   ```
   run_all.sh <build> <out> qmc_configs.txt qmc_density.txt   # ~5 min, 2e5 events per test
   python3 analyse.py <out>                                     # writes <out>/REPORT.md
   ```
   `qmc_configs.txt` / `qmc_density.txt` are plain-text conversions of
   `data/configurations/QMC_configs.out.gz` (header `N A`, then per configuration A lines
   `iso x y z` and a weight line; density: `r rho_p rho_n`).
2. Step 4(f), full kinematics against Achilles Veto/Continuous/Base, and Step 4(h) data runs:
   ```
   python3 compare_achilles.py <build> <out> Geant 40000
   python3 compare_achilles.py <build> <out> NN 60000
   python3 tabulate.py <out>/compare_NN.json NN
   python3 data_runs.py <build> <out> 100000 [alg,alg,...] [json name] [NN|Geant|NASA]
   python3 plots.py <forward out> <compare_NN.json> <data_runs.json> carlson.csv transparency.csv <figs>
   python3 plots_published.py <nasa.json> <nasa_sudakov.json> <geant.json> carlson.csv transparency.csv <figs>
   RANGES=transp python3 data_runs.py ...                       # one range only
   python3 chi2_table.py carlson.csv transparency.csv <runs.json> [...] [--json out.json]
   python3 plots_step2.py <val dir> <figs>                      # Pauli / coherent / in-medium / acceptance
   python3 plots_ipm.py <val dir> <density> <figs>              # IPM-corrected coherent veto
   python3 coherent_oracle.py <density>                         # IPM eikonal oracle, no Pauli
   ```
   Card keys in `compare_achilles.py`: `_cyl` sets `Probability: Cylinder`, `_med` sets
   `InMedium: NonRelativistic` (which the Sudakov tracks also use); in `data_runs.py`, `_pot`
   sets `PotentialProp: True`, `_cap1000` raises the cross-section cap, and `{ipm}` cards get the
   `ipm_veto.py` table for the run's density.
   Transparency mode also writes `<name>_accepted.txt`: events whose leading final-state nucleon of
   the kicked species lost < `accept_dE` (80 MeV) and changed momentum by < `accept_dp`
   (300 MeV/c); set both under `Cascade/Params`.

Runcard block (under `Cascade:`), all keys optional:
```
  Algorithm: SudakovB          # or SudakovA
  Sudakov:
    Law: Eikonal               # B: Eikonal (default) | Kappa1 (Achilles Gaussian)
    X0: 2.6                    # B: eikonal opacity
    Profile: Gauss             # B: Gauss | Box | Atom
    Headroom: 1.3              # majorant factor c >= 1 (needed for reweighting)
    Shells: 16                 # A: radial shells of the density majorant
    FoldWidth: 0.8             # A: 3D Gaussian folding of rho (fm)
    FiniteA: true              # A: (A-m)/A depletion and IPM void factor
    Pauli: Step                # Step (local Fermi gas) | FermiDirac | None (blocking off)
    PauliAttempt: Single       # Single (one attempt per pair, as Achilles; default) | Retry (hazard thinning)
    CoherentVeto: false        # true: veto Pauli-allowed elastic NN final states with prob. F_t(q)^2
    CoherentTable: {q: [...], v: [...]}  # instead of F^2: tabulated veto v(q) (q in MeV/c);
                               # python3 ipm_veto.py <density> prints the IPM table for 12C
    MediumMax: 1.0             # A: bound on the InMedium factor in the majorant
    PauliDelta: 0.05
    SigmaCap: 500              # mb; cap of the NN cross section (see AchillesCollision)
    Variations:                # weights for sigma -> s * sigma
      - {Name: sig08, SigmaScale: 0.8}
```
