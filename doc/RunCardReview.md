# Run Card, Parsing, and Rules — Review and Improvement Plan

This document reviews how Achilles sets up and parses the run card, how the
validation rules are implemented, and what works well versus what does not. It
ends with a concrete plan to better couple the `NuclearModels` and `Nuclei`
sections so that **multiple spectral functions can be used for different
nuclei in a single run**.

All file/line references are against the state of the tree when this review was
written.

---

## 1. How it works today

### 1.1 Run card structure

A run card (`run.yml`, `data/default/run.yml`) is a YAML document with these
top-level sections:

| Section          | Purpose                                                        |
|------------------|----------------------------------------------------------------|
| `Main`           | NEvents, cuts toggles, output format/name/zip, decays          |
| `SherpaOptions`  | BSM/Sherpa model + param card                                  |
| `Processes`      | list of leptonic in/out states                                |
| `Beams`          | list of beams (PID + beam params)                             |
| `Cascade`        | intranuclear cascade settings                                 |
| `NuclearModels`  | **list** of hard-interaction nuclear models                   |
| `Nuclei`         | **list** of target nuclei (each usually an `!include`)        |
| `HardCuts`       | generation-level cuts                                          |
| `Options`        | `Initialize` (seed/accuracy) + `Unweighting`                  |
| `Backend`        | cross-section backend                                          |

Nucleus definitions are normally pulled in with the custom `!include` tag,
e.g. `data/default/12C.yml`, `16O.yml`, `40Ar.yml`, `1H.yml`, `1N.yml`.

### 1.2 Parsing — the `Settings` class

`include/Achilles/Settings.hh`, `src/Achilles/Settings.cc`.

* `Settings(filename, rules="data/Rules.yml")` loads the YAML, recursively
  resolves `!include` tags (`IncludeFile`, Settings.cc:351), loads the rule
  file, and runs validation. `EventGen` holds a `Settings config;`
  (`EventGen.hh:56`) built from the run card, so **every run is validated
  against `data/Rules.yml`**.
* Access is via a path mini-language: `settings["Main/Output/Format"]`,
  `settings.GetAs<T>("Cascade/Step")`. `ParsePath` (Settings.cc:361) splits on
  `/`, turns all-digit segments into sequence indices, and supports a `*`
  wildcard. `seek` (Settings.cc:380) walks maps/sequences and can return
  multiple matches for a wildcard.
* `Exists` is implemented by catching the "not found" exception from
  `operator[]` (Settings.cc:345).

### 1.3 Rules / validation

`data/Rules.yml` (and the cascade-only `data/hA_Rules.yml`) drive
`SettingsValidator` (Settings.cc:46–294). Three rule types exist:

1. **`RequiredFields`** — each path in `Fields` must `Exists`
   (`ValidateRequiredFields`, Settings.cc:84).
2. **`ConditionalInteractionOption`** — if certain cascade interactions are
   present, a target interaction must carry a required option value
   (`ValidateConditionalInteractionOption`, Settings.cc:95).
3. **`RangeConstraint`** — min/max on a scalar `Path`, or on a `Parameter`
   inside each matching element of a `List` (Settings.cc:167–294).

### 1.4 Nuclear models

`src/Achilles/NuclearModel.cc`, `include/Achilles/NuclearModel.hh`.

* `LoadModels(const Settings&)` (NuclearModel.cc:193) iterates
  `settings["NuclearModels"]`, instantiates each via the
  `NuclearModelFactory`, and stores them in a
  `ModelMap = unordered_map<NuclearMode, unique_ptr<NuclearModel>>`
  **keyed by interaction mode** (Quasielastic, Resonance, MEC, …). Two models
  with the same `NuclearMode` is a fatal error (NuclearModel.cc:211).
* `QESpectral` reads its spectral functions directly from the model node:
  ```cpp
  spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
  spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()}
  ```
  (NuclearModel.cc:542). `FortranModel` instead reads a `ConfigFile`
  (`FNuclearModel.cc:34`) whose contents are the proton/neutron SF paths
  (e.g. `data/info_C12_pke.data` → `pke12{p,n}_tot.data`).
* `SpectralFunction` (`src/Achilles/SpectralFunction.cc`) just reads a numeric
  `(ne, np)` grid. **It carries no nucleus identity (no Z, A) at all** — the
  nucleus is implied only by which file path the user supplied.

### 1.5 Nuclei

* `EventGen` reads the whole list:
  `nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei")`
  (EventGen.cc:61).
* The decode lives in `convert<std::vector<std::shared_ptr<Nucleus>>>`
  (`Nucleus.hh:321`) and the per-nucleus `convert<Nucleus>` (`Nucleus.hh:288`),
  which special-cases `1H`/`1N` and otherwise builds densities, Fermi gas, and
  potential.

### 1.6 The model ↔ nucleus coupling — the core problem

`EventGen` constructor (EventGen.cc:101):

```cpp
for(const auto &nucleus : nuclei) {
    auto models = LoadModels(config);            // (A) same global config every time
    for(auto &model : models) {
        auto groups = ProcessGroup::ConstructGroups(config, model.second.get(),
                                                     beam, nucleus);
        for(auto &group : groups) {
            ...
            group.second.SetupBackend(config, std::move(model.second), p_sherpa); // (B)
            ...
        }
    }
}
```

Two things stand out:

* **(A) The model is built from the global `NuclearModels` section, not from
  anything tied to `nucleus`.** `ConstructGroups` (Process.cc:308) accepts both
  a `NuclearModel*` and a `Nucleus`, but the model's spectral function was
  already fixed by the run card path. The nucleus PID reaches the model only
  through `PhaseSpace(nuc_id)` (wired via Process.cc:357 →
  `XSecBackend::SetupChannels` → `Channels.cc` `model->PhaseSpace(nuc_id)`),
  which for `QESpectral` merely flips the `is_hydrogen` / `is_free_neutron`
  flags (NuclearModel.cc:645). **So if you list `12C`, `16O`, and `40Ar`
  together with a `QESpectral` model pointing at `pke12{p,n}_tot.data`, all
  three nuclei are generated with the carbon spectral function.** There is no
  error, just silently wrong physics. (An existing `// TODO` at
  NuclearModel.cc:475 — "make it work with multiple nuclei" — acknowledges
  this.)

* **Prior art: `Coherent` already binds to a nucleus.** The `Coherent` model
  stores an expected `nucleus_pid` from `config["NuclearModel"]["Nucleus"]` and
  `Coherent::PhaseSpace` *throws* `InvalidChannel` when the run's nucleus PID
  does not match (NuclearModel.cc:529). This is exactly the per-nucleus binding
  pattern the QE/RES spectral models lack — and it means a multi-nucleus run
  that includes a `Coherent` model breaks today unless every nucleus matches
  that one PID. The plan below generalizes Coherent's pattern rather than
  inventing a new one. Nucleus PIDs are the PDG nuclear codes
  `10LZZZAAAI` built in `Nucleus::Initialize` (e.g. `12C` → `1000060120`),
  so a model can unambiguously identify its intended target by Z and A.

* **(B) `std::move(model.second)` inside the inner `groups` loop.** `groups` is
  a `map` keyed by final-state multiplicity; if a single model yields more than
  one multiplicity group, the model `unique_ptr` is moved on the first group
  and the remaining groups receive a moved-from null pointer. This is a latent
  bug that only hides because today most QE/coherent models produce a single
  multiplicity.

---

## 2. The Good

* **Data-driven validation.** Rules live in YAML, not C++. Adding a required
  field or a numeric bound is a one-line edit, and the three rule primitives
  (`RequiredFields`, `ConditionalInteractionOption`, `RangeConstraint`) cover a
  lot of ground. The list-matching `RangeConstraint` (match-by-name, then check
  a sub-parameter) is genuinely nice.
* **`!include` composition.** Nucleus, Sherpa, interaction, and option blocks
  are factored into reusable files under `data/default/`. The recursive visitor
  resolves them transparently.
* **Path mini-language with wildcards.** `settings["A/B/2/C"]` and `*` wildcards
  give concise access and keep call sites readable.
* **List-based sections.** `Processes`, `Beams`, `NuclearModels`, and `Nuclei`
  are already sequences, so the schema is forward-compatible with multiplicity
  even where the runtime is not yet.
* **Helpful failure modes for some errors.** Unknown model/unweighter names
  print a "did you mean …" suggestion (NuclearModel.cc:220, Process.cc:329).

---

## 3. The Bad / Issues

### 3.1 Correctness / footguns

1. **No nucleus ↔ model / spectral-function binding** (Section 1.6). The single
   biggest issue, and the subject of the plan below. Multiple nuclei silently
   share one spectral function.
2. **`std::move` in a loop** (EventGen.cc:110) can null a model for
   multi-multiplicity groups.
3. **Spectral functions carry no identity.** `SpectralFunction` never records
   Z/A, so nothing can cross-check that `pke12*` was loaded for a carbon
   target. A typo in a path produces wrong physics, not an error. Worse, the SF
   filenames encode only the mass number `A` (`pke12*`, `pke40*`), not `Z`, so
   two isobars (same `A`, different `Z`) would silently share a file.
4. **`ConfigFile` is dead config for `QESpectral`.** The C++ `QESpectral`
   ignores `ConfigFile` entirely (only `SpectralP`/`SpectralN` are read), yet
   `run.yml` and `data/default/run.yml` both set `ConfigFile` on the
   `QESpectral` block. Only `FortranModel` reads it. This is confusing and
   invites copy-paste errors.

### 3.2 Schema drift / inconsistency

5. **Examples contradict the rules.** `data/Rules.yml` requires
   `Options/Unweighting` and `Options/Initialize`, but every card in
   `examples/` puts `Initialize:` and `Unweighting:` at the **top level** (e.g.
   `examples/run_electron_12C_SF_1000.yml:17`). Those cards would fail
   validation against the current rules — they are stale.
6. **README is out of date.** `README.md:135` documents singular `Process`,
   `Initialization`, `Unweighting`, `Nuclear Model`, and `Nucleus` sections and
   states "only a single isotope and nucleus is supported", while the live
   schema uses plural list sections (`Processes`, `NuclearModels`, `Nuclei`)
   and an `Options` block. `run_hnl.yml` is also on the old singular schema.
7. **Two parallel rule files** (`data/Rules.yml`, `data/hA_Rules.yml`) with the
   same `ConditionalInteractionOption`/`RangeConstraint` blocks duplicated
   verbatim. They will drift.

### 3.3 Validator robustness

8. **`spdlog::error(rule.error_message, …)` uses runtime strings as format
   strings.** `ValidateRequiredFields` (Settings.cc:88) passes the YAML
   `ErrorMessage` straight to fmt as a format string. With modern fmt this is a
   compile-time-checked API being fed a runtime string; a stray `{`/`}` in a
   message, or an arg-count mismatch, throws `fmt::format_error` instead of
   reporting the real problem. Use `fmt::runtime(...)` and validate arg counts.
9. **`ConditionalInteractionOption` only runs when `Cascade/Mode` is defined**
   (Settings.cc:99), but the standard cards use `Cascade/Run` +
   `Cascade/Interactions` and never set `Cascade/Mode`. So in practice this rule
   short-circuits to "valid" for the very cards that ship — the validation it
   promises is mostly not happening.
10. **Required cascade fields are not actually required.**
    `RangeConstraint` on `Cascade/Step` has `Optional: False` and a hard
    `min`, but if `Cascade/Run: False` the step is irrelevant; conversely there
    is no rule asserting `Cascade/Interactions` exists when `Cascade/Run:
    True`. Several real preconditions live in C++ exceptions, not rules.
11. **No "unknown key" detection.** A misspelled `NEvenst:` or `Sectral:` is
    silently ignored; defaults quietly take over. There is no schema closure.
12. **Numeric-only path segments are always indices.** `ParsePath`
    (Settings.cc:370) turns any all-digit key into a sequence index, so a YAML
    **map** key that happens to be numeric can never be addressed. Minor, but a
    sharp edge.

### 3.4 Usability

13. **Heavy redundancy per nucleus.** Each nucleus file repeats binding, Fermi
    momentum, density files, Fermi-gas block, and potential — most of which are
    identical boilerplate. `1N.yml` even has a stray `Mode: 3s` typo
    (`data/default/1N.yml:18`) that is silently ignored.
14. **Model/SF/nucleus must be kept in sync by hand across three places** (the
    `NuclearModels` SF path, the `ConfigFile`, and the chosen `Nuclei` entry)
    with no tooling to catch a mismatch.

---

## 4. Suggested improvements (streamlining)

Ordered roughly by value/effort.

1. **Fix the schema drift first** (low effort, high trust):
   * Regenerate `examples/*.yml` to put `Initialize`/`Unweighting` under
     `Options`, or relax the rule to accept either — pick one and make the
     repo consistent.
   * Update `README.md` and `run_hnl.yml` to the plural/`Options` schema and
     drop the "single nucleus only" claim once Section 5 lands.
   * Merge `Rules.yml` and `hA_Rules.yml` via `!include` of a shared block so
     the common rules exist once.

2. **Make the validator robust** (low effort):
   * Wrap runtime messages in `fmt::runtime(...)` and guard arg counts
     (Settings.cc:88, and the other `spdlog::error(rule.error_message, …)`
     sites).
   * Make `ConditionalInteractionOption` trigger on `Cascade/Run: True`
     (the field the cards actually use) rather than `Cascade/Mode`.
   * Add an **`OneOf` / `MutuallyExclusive`** rule type and a
     **`RequiredIf`** rule type (e.g. "`Cascade/Interactions` required if
     `Cascade/Run` is true"). These cover the preconditions currently buried in
     C++ throws.

3. **Add unknown-key warnings.** After load, walk the tree and warn on keys not
   present in a known-keys manifest (can start as warn-only). Catches the
   `NEvenst` class of bugs that defaults currently mask.

4. **Remove dead/confusing config.** Drop `ConfigFile` from `QESpectral`
   blocks (or make `QESpectral` honor it the way `FortranModel` does, and
   document one canonical way). Fix the `1N.yml` `Mode: 3s` typo.

5. **Reduce nucleus boilerplate.** Provide a `data/default/nuclei/` library of
   canonical nucleus files (already mostly there) plus a thin
   `NucleusDefaults.yml` that individual cards override, so a card only states
   what differs.

6. **Surface a `--validate`/dry-run flag** that loads + validates the card and
   prints the fully-resolved settings (`Settings::Print` already exists) without
   running, so users can debug cards cheaply.

---

## 5. Plan: tie `NuclearModels` to `Nuclei`, with per-nucleus spectral functions

**Goal.** Allow a single run to contain several nuclei, each with its own
spectral function(s) / model parameters, while keeping the common case (one
nucleus) as terse as today.

### 5.1 Target run-card schema

Introduce an explicit, optional binding between a nucleus and the models that
apply to it. Two complementary mechanisms:

**(a) Per-nucleus model/SF overrides (recommended primary mechanism).**

```yaml
NuclearModels:
  - NuclearModel:
      Model: QESpectral
      FormFactorFile: "FormFactors.yml"
      Ward: None
      # Default SF used if a nucleus does not override it:
      SpectralP: data/Spectral_Functions/pke12p_tot.data
      SpectralN: data/Spectral_Functions/pke12n_tot.data

Nuclei:
  - Nucleus:
      !include "data/default/12C.yml"
      # inherits the default carbon SF above
  - Nucleus:
      !include "data/default/40Ar.yml"
      Models:                      # NEW: per-nucleus overrides keyed by model
        QESpectral:
          SpectralP: data/Spectral_Functions/pke40p_tot.data
          SpectralN: data/Spectral_Functions/pke40n_tot.data
```

**(b) Explicit selection by name (for advanced multi-model setups).** Give each
model an optional `Label`, and let a nucleus list which labels it uses:

```yaml
NuclearModels:
  - NuclearModel: { Label: QE_C12, Model: QESpectral, SpectralP: .../pke12p..., ... }
  - NuclearModel: { Label: QE_Ar40, Model: QESpectral, SpectralP: .../pke40p..., ... }

Nuclei:
  - Nucleus: { !include "data/default/12C.yml",  Models: [QE_C12] }
  - Nucleus: { !include "data/default/40Ar.yml", Models: [QE_Ar40] }
```

Mechanism (a) keeps the single-nucleus card unchanged (full backward
compatibility); mechanism (b) lifts the current "one model per `NuclearMode`"
restriction by scoping that restriction **per nucleus** instead of globally.

### 5.2 Code changes

1. **`LoadModels` becomes nucleus-aware.** Change the free function to
   `LoadModels(const Settings&, const Nucleus&)` (or
   `LoadModelsFor(nucleus_node, global_models_node)`), so it can merge the
   global model node with the nucleus's `Models:` override before constructing.
   Concretely: deep-copy the model's YAML node, overlay the per-nucleus keys,
   then call the factory. Keep the existing global path as the default branch.
   *Files:* `NuclearModel.cc:193`, `NuclearModel.hh:151`.

2. **Scope the duplicate-mode check per nucleus.** The "Multiple nuclear models
   for mode" error (NuclearModel.cc:211) should be evaluated against the model
   set selected for *that* nucleus, not globally, so two `QESpectral` configs
   for different nuclei stop colliding.

3. **Rework the `EventGen` loop** (EventGen.cc:101). Build the model map
   *inside* the per-nucleus loop from the merged config, and **fix the
   `std::move` bug**: move the model into the backend exactly once, or clone per
   group. Sketch:
   ```cpp
   for(const auto &nucleus : nuclei) {
       auto models = LoadModels(config, *nucleus);   // merged/overridden
       for(auto &[mode, model] : models) {
           auto groups = ProcessGroup::ConstructGroups(config, model.get(), beam, nucleus);
           // hand each group its own model instance (move once / clone)
           for(auto &[mult, group] : groups) {
               group.SetupBackend(config, /* one model per group */, p_sherpa);
               process_groups.push_back(std::move(group));
           }
       }
   }
   ```
   This needs a way to give every group an owning model — either construct one
   model per (nucleus, mode, multiplicity) up front, or add a `Clone()` to
   `NuclearModel`.

4. **Give `SpectralFunction` an identity and validate it.** Add optional
   `(Z, A)` metadata to `SpectralFunction` (passed in at construction, or read
   from a header line in the SF file). At model build time, assert the SF's
   `(Z, A)` matches the owning `Nucleus`. Emit a clear error on mismatch instead
   of silently generating wrong physics. *Files:* `SpectralFunction.hh/.cc`,
   model constructors (NuclearModel.cc:542, :659; `FNuclearModel.cc`).

5. **New validation rules** in `data/Rules.yml`:
   * `RequiredFields` per QE/RES spectral model: each `QESpectral` block (or its
     per-nucleus override) resolves to readable `SpectralP`/`SpectralN`.
   * A `CrossReference` rule type that checks every `Models:` label referenced
     by a nucleus exists in `NuclearModels` (mechanism b).
   * Optional `RequiredIf`: if a nucleus is not `1H`/`1N`, a QE model must
     provide both proton and neutron SFs.

6. **Fortran models.** `FortranModel` already takes a `ConfigFile` that lists
   the SF paths, so per-nucleus support is just allowing the nucleus override to
   swap `ConfigFile` (e.g. `info_C12_pke.data` vs `info_Ar40_pke.data`). Make
   sure the merged node flows into `FNuclearModel.cc:34`.

### 5.3 Migration / compatibility

* The default branch of `LoadModels` (no `Models:` override, no `Label`) must
  reproduce today's behavior byte-for-byte, so existing single-nucleus cards
  and the `examples/` set keep working.
* Add at least one multi-nucleus example
  (e.g. `examples/run_multinuc_C12_Ar40.yml`) and a regression test analogous to
  `test/test_nuclear_model.cc` asserting that each nucleus picks up its own SF
  (e.g. by checking the SF normalization or a sampled `S(p,E)` value differs
  between the C and Ar groups).
* Update `README.md` §"Run card" to document the new `Models:` override and drop
  the single-nucleus limitation.

### 5.4 Suggested sequencing

1. Schema-drift cleanup + validator-robustness fixes (Section 4.1–4.2) — safe,
   immediately useful, unblocks trustworthy testing.
2. `SpectralFunction` identity + mismatch error (4 above) — turns the silent
   footgun into a hard error even before multi-nucleus lands.
3. Nucleus-aware `LoadModels` + `EventGen` loop fix (1–3 above) with the
   per-nucleus override schema (mechanism a).
4. Named-model selection (mechanism b) + `CrossReference` rule, if multi-model
   workflows need it.
5. Docs + multi-nucleus example + regression test.
