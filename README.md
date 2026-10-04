# NitroAfrica modifications to eCLM

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21490379.svg)](https://doi.org/10.5281/zenodo.21490379)

This repository contains the modified eCLM source code developed for the
NitroAfrica study of atmospheric nitrogen deposition composition and soil
nitrogen emissions in African savanna ecosystems.

# eCLM

[![CI](https://github.com/HPSCTerrSys/eCLM/actions/workflows/CI.yml/badge.svg)](https://github.com/HPSCTerrSys/eCLM/actions/workflows/CI.yml)
[![docs](https://github.com/HPSCTerrSys/eCLM/actions/workflows/docs.yml/badge.svg)](https://github.com/HPSCTerrSys/eCLM/actions/workflows/docs.yml)
[![latest tag](https://badgen.net/github/tag/HPSCTerrSys/eCLM)](https://github.com/HPSCTerrSys/eCLM/tags)

eCLM is based from [Community Land Model 5.0 (CLM5.0)]. It has the same modelling capabilities as CLM5 but with a more simplified infrastructure for build and namelist generation. The build system is handled entirely by CMake and namelists are generated through a small set of Python scripts. Only Fortran source codes necessary for a functional land model simulation were imported from CLM5. 

Unlike CLM5, there are no built-in batch scripts in eCLM. It is up to system maintainers or users to craft their own workflows by combining the basic tools in this repo plus the native tools in their respective platforms.

> [!WARNING]
> eCLM is still experimental and being extensively tested. Use it at your own risk!

## Reproducible release

The source-code release archived for peer review is
[`v0.1.3-peer-review`](https://github.com/profdam/eCLM-NitroAfrica/tree/v0.1.3-peer-review).
It corresponds to commit
`61e32b171305694767737a5a84696f335c38cf01` and the version-specific Zenodo
DOI [10.5281/zenodo.21490379](https://doi.org/10.5281/zenodo.21490379).

```bash
git clone https://github.com/profdam/eCLM-NitroAfrica.git
cd eCLM-NitroAfrica
git checkout v0.1.3-peer-review
git rev-parse HEAD
```

The final command should return the commit recorded above. DOI
`10.5281/zenodo.21490378` is the concept DOI for all versions. Use the
version-specific DOI when reproducing or citing this release.

## Build

eCLM uses CMake and requires compatible C and Fortran compilers, MPI, NetCDF,
HDF5 and PnetCDF. From the repository root, configure an out-of-source build:

```bash
export ECLM_INSTALL="$PWD/install"
cmake -S src -B build \
  -DCMAKE_INSTALL_PREFIX="$ECLM_INSTALL" \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpifort
cmake --build build
cmake --install build
```

Compiler and dependency discovery is platform-specific. Consult the
[bundled build guide](docs/users_guide/building_eCLM/README.md) and the
[upstream eCLM documentation](https://hpscterrsys.github.io/eCLM) for the
available compiler and platform options. The exact compiler versions, library
versions and CMake options used for the NitroAfrica production simulations are
not recorded in this repository and should not be inferred from the generic
examples.

The namelist generator requires Python 3 and can be installed from the clone:

```bash
python3 -m pip install --user ./namelist_generator
```

## NitroAfrica deposition input

The modified deposition stream requires chemically resolved nitrogen
deposition. Its default variable names are:

- `NDEP_NHx_month` for reduced nitrogen deposition;
- `NDEP_NOy_month` for oxidised nitrogen deposition.

NHx is routed to the soil mineral NH4+ pool and NOy to the soil mineral NO3-
pool. Before a production run, check the units and temporal interpretation of
both variables against `src/clm5/main/ndepStreamMod.F90`. Further scientific
and technical details are provided in [MODIFICATIONS.md](MODIFICATIONS.md).

## Inputs and configuration needed for reproduction

The repository contains the modified source code, but not all files required
to reproduce the Dahra, Lamto or Korhogo experiments. Users also require the
applicable domain and surface datasets, atmospheric and deposition forcing,
initial or restart conditions, namelists and history-field configuration.
Large model outputs and restricted observations are not included.

The repository does not currently provide persistent identifiers or checksums
for all of those inputs. Consequently, the tagged source code is citable and
archived, but the complete manuscript workflow cannot be reconstructed from
this repository alone.

## Minimum verification

Before running a full experiment:

1. build from a clean directory;
2. run a short test simulation;
3. confirm that `NDEP_NHx_month` and `NDEP_NOy_month` are read correctly;
4. confirm that NHx is routed to NH4+ and NOy to NO3-;
5. inspect the nitrogen-balance diagnostics;
6. verify native output units and time integration; and
7. record the Git commit, compiler and dependency versions.

## Citation

Please cite the version-specific archived release:

> Dahunsi, A. M., Delon, C., & Solmon, F. (2026). *NitroAfrica modifications
> to eCLM* (Version v0.1.3-peer-review) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.21490379

[Community Land Model 5.0 (CLM5.0)]: https://github.com/ESCOMP/CTSM/tree/release-clm5.0
