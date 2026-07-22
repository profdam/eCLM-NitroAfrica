# NitroAfrica modifications to eCLM

## Overview

This repository contains a modified version of the eCLM land-surface model developed for the NitroAfrica project. The modifications were introduced to represent the chemical composition of atmospheric nitrogen deposition and to improve the simulation and diagnosis of soil nitrogen emissions in African savanna ecosystems.

The development focuses on the Dahra (Senegal), Lamto and Korhogo (Cote d'Ivoire) study sites and includes changes related to ammonium and oxidised-nitrogen deposition, soil nitric oxide emissions, nitrous oxide diagnostics, ammonia volatilisation, rainfall-induced emission pulses and canopy reduction of soil emissions.

The repository contains the modified source code required for the simulations reported in the associated manuscript. Large forcing datasets, surface datasets, model output, restart files and compiled executables are not included.

## Upstream model

The modified model is derived from the following upstream eCLM version:

* Model: eCLM
* Upstream repository: `HPSCTerrSys/eCLM`
* Upstream branch: `master`
* Base commit: `9dd644fbd0c532fdf33cdab33f17cfa2160b9fbf`
* Date retrieved: 3 March 2025

The original eCLM source-code history, copyright statements and licence conditions are retained.

## Purpose of the modifications

The original model represents total atmospheric nitrogen deposition as a single nitrogen flux. The modifications introduced here distinguish between reduced and oxidised forms of deposited nitrogen and route them to the corresponding soil mineral nitrogen pools.

The model was also extended to diagnose soil emissions of NH3 and NO using process formulations adapted for studies of nitrogen cycling in tropical and semi-arid ecosystems.

The principal objectives of the modifications are to:

1. distinguish NHx and NOy atmospheric deposition;
2. route deposited NHx to the soil NH4+ pool and deposited NOy to the soil NO3− pool;
3. retain total nitrogen deposition as the sum of the NHx and NOy components;
4. add alternative process-based diagnostics for soil NO, N2O and NH3 emissions;
5. represent rainfall-induced pulses in soil NO emissions;
6. account for canopy reduction of soil NO emissions;
7. distinguish potential soil emissions from emissions reaching the atmosphere;
8. expose the additional nitrogen fluxes and state variables in model history output.

## Scientific modifications

### 1. Separation of NHx and NOy deposition

The atmospheric nitrogen-deposition input was separated into two components:

* `NDEP_NHx_month`, representing reduced nitrogen deposition;
* `NDEP_NOy_month`, representing oxidised nitrogen deposition.

The model retains the total deposition flux as:

Total nitrogen deposition = NHx deposition + NOy deposition

Separate forcing variables were added to the atmosphere-to-land data structure:

* `forc_nhxdep_grc`;
* `forc_noydep_grc`.

The existing total-deposition variable, `forc_ndep_grc`, is calculated as the sum of these two components.

### 2. Deposition input-stream handling

The nitrogen-deposition stream reader was modified to read two deposition variables from the same input stream.

The unit-checking procedure was extended to:

* inspect the units of the NHx and NOy variables independently;
* apply the required temporal conversion to each component;
* retain support for deposition data supplied either as annual totals or as fluxes;
* calculate total deposition after unit conversion.

When deposition is supplied through the coupler, the NHx and NOy fluxes are also transferred separately and converted from kg N m−2 s−1 to g N m−2 s−1 before being passed to the land model.

### 3. Routing of deposited nitrogen to soil pools

The deposition fluxes were routed according to their chemical form:

* deposited NHx is added to the soil mineral NH4+ pool;
* deposited NOy is added to the soil mineral NO3− pool.

The following column-level flux variables were introduced:

* `nhxdep_to_sminnh4_col`;
* `noydep_to_sminno3_col`.

The original total deposition variable, `ndep_to_sminn_col`, is retained for nitrogen-balance calculations and diagnostic purposes.

This treatment replaces the original assumption that all deposited nitrogen enters the soil NH4+ pool.

### 4. Soil NO and NOx emissions

Additional calculations were introduced to represent soil NOx production associated with nitrification and denitrification.

The implementation distinguishes:

* NOx produced by nitrification;
* NOx produced by denitrification;
* potential NOx production;
* NOx remaining after process constraints;
* NOx reaching the atmosphere after canopy reduction.

Vertically resolved and column-integrated fluxes are retained separately.

The NOx:N2O partitioning follows the alternative formulation implemented for this study, based on soil environmental controls and the relative contributions of nitrification and denitrification.

### 5. Rainfall-induced NO emission pulses

A rainfall-pulse treatment was introduced for soil NO production.

The implementation tracks antecedent soil dryness using vertically resolved variables representing:

* the duration of the dry period;
* the rainfall-pulse enhancement factor.

The pulse formulation uses:

* a volumetric soil-moisture threshold of 0.175 m3 m−3;
* antecedent dry-period duration;
* a rainfall-triggered enhancement factor;
* temporal decay of the pulse after wetting.

The pulse is applied to the realised soil NOx production after nitrogen competition has been resolved.

The following variables were added:

* `ldry_vr_col`;
* `pfactor_vr_col`;
* `pulse_fac_col`.

### 6. Canopy reduction of soil NO emissions

A canopy reduction factor was introduced to distinguish soil NO production from the fraction of the flux that reaches the atmosphere.

The reduction factor is calculated from exposed leaf and stem area indices:

CRF = 0.5 × [exp(−ks × SAI) + exp(−kc × LAI)]

where:

* `LAI` is exposed leaf area index;
* `SAI` is exposed stem area index;
* `kc = 0.32`;
* `ks = 11.6`.

The canopy reduction factor is constrained to the interval from zero to one.

Patch-level canopy reduction factors are mapped to model columns before being applied to the soil NO fluxes.

The implementation therefore distinguishes between:

* soil NO or NOx production below the canopy;
* NO or NOx emissions reaching the atmosphere.

### 7. Alternative N2O diagnostics

Additional N2O flux variables were introduced to diagnose emissions associated separately with:

* nitrification;
* denitrification;
* total soil N2O production.

The alternative diagnostics are retained alongside the original CLM5 N2O calculation. The original CLM5 N2O flux is therefore not replaced.

The additional variables permit direct comparison between the default CLM5 formulation and the alternative process representation used in the NitroAfrica analysis.

### 8. NH3 volatilisation

A soil NH3 volatilisation treatment was incorporated into the soil-biogeochemistry calculations.

The implementation includes:

* potential NH3 volatilisation;
* realised NH3 volatilisation;
* vertically resolved NH3 fluxes;
* column-integrated NH3 fluxes;
* soil aqueous or available NH3 diagnostics;
* partitioning between atmospheric emission and canopy capture;
* dependence on soil temperature;
* dependence on soil NH4+ availability;
* dependence on soil pH;
* adsorption effects related to soil properties;
* aerodynamic and friction-velocity controls.

The following categories of variables were introduced:

* potential NH3 volatilisation;
* actual NH3 volatilisation;
* NH3 emitted to the atmosphere;
* NH3 intercepted by the canopy;
* soil NH3 state diagnostics.

For the experiments reported in the associated study, the relevant soil-pH assumptions and parameter values are documented in the manuscript and model configuration.

### 9. Potential and realised nitrogen fluxes

Several nitrogen-emission variables are retained in both potential and realised forms.

Potential fluxes represent the unconstrained process rate before all competition, availability, canopy and environmental limitations are applied.

Realised fluxes represent the final modelled flux following the relevant constraints.

This distinction was introduced for:

* NH3 volatilisation;
* NOx production from nitrification;
* NOx production from denitrification;
* above-canopy NOx emissions.

### 10. Nitrogen state and flux data structures

The nitrogen-state and nitrogen-flux data structures were extended to allocate, initialise, update and archive the new variables.

The additions include:

* separate NHx and NOy deposition fluxes;
* vertically resolved and integrated NH3 fluxes;
* vertically resolved and integrated NOx fluxes;
* nitrification and denitrification components;
* potential and realised emission fluxes;
* above-canopy emission fluxes;
* rain-pulse state variables;
* diagnostic NO:N2O partitioning variables;
* total soil NOx and N2O diagnostics;
* crop-specific nitrogen-emission diagnostics where applicable.

### 11. History output

The model history-output registration was extended to expose the additional nitrogen variables for analysis.

The new or revised history fields include variables representing:

* NHx deposition;
* NOy deposition;
* total nitrogen deposition;
* NH3 volatilisation;
* potential NH3 volatilisation;
* NH3 emitted to the atmosphere;
* NH3 intercepted by the canopy;
* soil NH3;
* NOx production from nitrification;
* NOx production from denitrification;
* potential NOx production;
* above-canopy NOx emissions;
* nitrification-derived N2O;
* denitrification-derived N2O;
* total soil NOx;
* total soil N2O;
* rainfall-pulse factors;
* canopy reduction factors;
* crop-specific nitrogen fluxes where activated.

The exact variable names available in an output file depend on the requested history-field list and model configuration.

## Modified source files

The following tracked source files contain the modifications associated with this version.

### Biogeochemistry driver and deposition routing

* `src/clm5/biogeochem/CNBalanceCheckMod.F90`
* `src/clm5/biogeochem/CNDriverMod.F90`
* `src/clm5/biogeochem/CNNDynamicsMod.F90`
* `src/clm5/biogeochem/CNVegetationFacade.F90`

### Atmosphere–land coupling and deposition input

* `src/clm5/cpl/lnd_import_export.F90`
* `src/clm5/main/atm2lndType.F90`
* `src/clm5/main/clm_driver.F90`
* `src/clm5/main/ndepStreamMod.F90`

### Soil nitrogen transformations and emissions

* `src/clm5/soilbiogeochem/SoilBiogeochemCompetitionMod.F90`
* `src/clm5/soilbiogeochem/SoilBiogeochemNStateUpdate1Mod.F90`
* `src/clm5/soilbiogeochem/SoilBiogeochemNitrifDenitrifMod.F90`
* `src/clm5/soilbiogeochem/SoilBiogeochemNitrogenFluxType.F90`
* `src/clm5/soilbiogeochem/SoilBiogeochemNitrogenStateType.F90`


## Study configurations

The modified model was developed for simulations at the following African savanna sites:

* Dahra, Senegal;
* Lamto, Côte d’Ivoire;
* Korhogo, Côte d’Ivoire.

The model is also being used in the development of a regional West African configuration.

The source-code repository does not automatically contain all namelists, forcing files, surface datasets, domain files or deposition datasets required to reproduce every experiment. These materials should be provided separately where redistribution is permitted.

## Input data

The modified source code requires deposition input containing separate NHx and NOy variables.

The expected default variable names in the modified deposition stream are:

NDEP_NHx_month
NDEP_NOy_month


Users applying different variable names must modify the relevant stream configuration or source-code defaults.

The unit attributes of both variables must be consistent with the unit-handling options supported by `ndepStreamMod.F90`.

## Compilation and testing

This modified version has been compiled and executed within the eCLM build system used for the NitroAfrica simulations.

Before using the code on another computing platform, users should:

1. configure the compiler, MPI, NetCDF, HDF5 and PnetCDF dependencies required by eCLM;
2. rebuild the model from a clean build directory;
3. verify that the modified deposition fields are read correctly;
4. confirm that NHx and NOy are routed to the intended mineral nitrogen pools;
5. perform a short simulation before starting a full experiment;
6. inspect nitrogen-balance diagnostics;
7. verify the units and temporal integration of all new output fluxes.

## Optional debugging diagnostics

The source code contains some commented diagnostic statements used during development to verify:

* NHx and NOy input-stream values;
* coupler transfer of deposition components;
* routing of deposited nitrogen to NH4+ and NO3−;
* nitrogen mass-balance terms.

These statements are inactive in the released configuration. They may be uncommented temporarily for debugging but should not be enabled during production simulations because they can produce large log files and increase input/output overhead.

## Repository scope

This repository is intended to contain:

* the original eCLM source code;
* the documented NitroAfrica source-code modifications;
* the original eCLM licence;
* this modification record;
* selected model configuration files where redistribution is permitted;
* scripts required to build or run the manuscript version;
* post-processing scripts required to reproduce reported diagnostics.

The repository is not intended to contain:

* model history files;
* restart files;
* compiled executables;
* object or module files;
* static or shared libraries;
* large NetCDF forcing files;
* restricted observational datasets;
* machine-specific temporary files.

## Known limitations

1. The added NO, N2O and NH3 formulations are diagnostic or experimental extensions developed for the NitroAfrica application.

2. The parameter values used in the emission formulations may require recalibration before application to ecosystems substantially different from the study sites.

3. The canopy reduction treatment represents bulk attenuation and does not explicitly resolve within-canopy chemistry or turbulent transport.

4. The rainfall-pulse treatment depends on model soil moisture and antecedent dry-period calculations. Its behaviour is sensitive to model hydrology and temporal resolution.

5. Soil pH is prescribed or externally constrained according to the model configuration. Dynamic soil-pH processes are not fully represented in the NitroAfrica experiments.

6. The distinction between NHx and NOy deposition depends on the availability and quality of chemically resolved deposition input data.

7. The modified code should not be assumed to reproduce the default eCLM solution when the new processes or deposition inputs are active.

8. Model outputs should be checked carefully for their native units before temporal integration or conversion to kg N ha−1 yr−1.

## Contributors

The modifications documented in this repository were developed by:

* Adeola Michael Dahunsi;
* Claire Delon;
* Fabien Solmon.

## Citation

Users of this modified version should cite:

1. the relevant CLM5 model description;
2. the original eCLM model and repository;
3. the NitroAfrica manuscript describing these modifications;
4. the archived software release and DOI, once available;
5. the scientific publications underlying the implemented NO, N2O, NH3, rainfall-pulse and canopy-reduction formulations.

The final manuscript citation and repository DOI will be inserted here after publication.


## Licence

This repository is derived from `HPSCTerrSys/eCLM`.

The original copyright notices, licence conditions and disclaimer are retained in the repository’s `LICENSE` file.

Redistribution and modification of the source code must comply with that licence. Source redistributions must retain the original copyright notice, licence conditions and disclaimer.

The names of the original copyright holders and contributors must not be used to endorse or promote derived products without prior written permission.

Any third-party datasets, parameter files or software components distributed separately remain subject to their respective licences.

## Disclaimer

This software is provided for research purposes without warranty. Users are responsible for checking the scientific suitability of the modified formulations, validating the model for their application and confirming compliance with the licences of all model inputs and dependencies.
