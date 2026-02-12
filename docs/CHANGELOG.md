# Changelog

All notable changes to the SIMC C++/ROOT project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- HMS spectrometer implementation (Phase 5c.4)
- Complete event generation for all reaction types
- Two-arm coincidence timing
- Full validation suite

---

## [0.5.3] - 2026-02-11

### Added - Phase 5c.3: Acceptance-Based Event Generation
- Acceptance-based angle generation for elastic H(e,e'p)
- Generate hadron angles within spectrometer acceptance (±60 mrad)
- Physics angle calculation from spectrometer angles
- Integration with SHMS transport

### Fixed
- Infinite loop bug: Removed overly tight kinematic constraint
- Event generation now achieves 100% efficiency
- Proper handling of spectrometer/physics angle conversion

### Validated
- 10,000 events: 100% generation efficiency
- SHMS acceptance: 17.1% (correct for config mismatch)
- Matches Fortran SIMC behavior

---

## [0.5.2] - 2026-02-11

### Added - Phase 5c.2: SHMS Spectrometer (Complete)
- Full SHMS class implementation
- 32 aperture planes:
  - HB (Heavy gas Cherenkov): 4 planes
  - Q1 (Quadrupole 1): 5 planes
  - Q2 (Quadrupole 2): 5 planes
  - Q3 (Quadrupole 3): 5 planes
  - Dipole: 13 planes (entrance, flare, 7 internal, exit)
- 5th order polynomial forward matrix
- Reconstruction matrix
- TRANSPORT coordinate system support
- Matrix file parsing (`shms_forward.dat`, `shms_recon.dat`)

### Changed
- Build system: Added `simc_shms` library
- CMakeLists.txt: Separated SHMS into own module
- Main loop: Integrated SHMS transport

### Validated
- Standalone test: 72.5% acceptance (10,000 events)
- Integrated test: 17.1% acceptance (explained by momentum mismatch)
- Transport algorithm matches Fortran SIMC

---

## [0.5.1] - 2026-02-XX

### Added - Phase 5c.1: Core Monte Carlo Transport
- `MonteCarloTransport` class
- Basic track state propagation
- Energy loss corrections in target
- Coordinate transformation utilities
- Placeholder acceptance checking

---

## [0.5.0] - 2026-02-XX

### Added - Phase 5b: Multiple Scattering
- Lynch-Dahl multiple scattering algorithm
- Beam multiple scattering through target
- Particle multiple scattering (electron, hadron)
- Material-dependent scattering calculations
- Integration with target geometry

### Files
- `include/MonteCarloTransport.h`
- `src/simc_core/MonteCarloTransport.cpp`

---

## [0.4.1] - 2026-02-XX

### Added - Phase 5a: Cross Sections
- `CrossSection` class
- Elastic H(e,e'p) cross section (dipole form factors)
- Quasi-elastic D(e,e'p) cross section
- Mott cross section calculation
- Form factor parameterizations

### Files
- `include/CrossSection.h`
- `src/simc_core/CrossSection.cpp`

---

## [0.4.0] - 2026-01-XX

### Added - Phase 4: Radiative Corrections
- `RadiativeCorrections` class
- Mo & Tsai formalism implementation
- Internal bremsstrahlung
- External bremsstrahlung
- Schwinger correction
- Peaking approximation
- Radiative tail generation

### Files
- `include/RadiativeCorrections.h`
- `src/simc_core/RadiativeCorrections.cpp`

---

## [0.3.0] - 2025-12-XX

### Added - Phase 3: Configuration System
- `ConfigManager` class
- JSON configuration file support
- Type-safe parameter access
- Default value handling
- Validation at load time
- Replaced Fortran namelists with modern JSON

### Files
- `include/ConfigManager.h`
- `src/simc_core/ConfigManager.cpp`
- `data/config/default.json`

### Dependencies
- Added nlohmann/json (header-only library)

---

## [0.2.0] - 2025-11-XX

### Added - Phase 2: Data Structures
- `SimcEvent` class (ROOT-compatible event structure)
- `TargetInfo` structure
- `SpectrometerInfo` structure
- `MainEvent` structure
- Physics constants (`PhysicsConstants.h`)
- ROOT dictionary generation for I/O

### Files
- `include/SimcEvent.h`
- `include/SimcTypes.h`
- `include/PhysicsConstants.h`

---

## [0.1.0] - 2025-10-XX

### Added - Phase 1: Project Setup
- CMake build system (modern target-based)
- ROOT framework integration
- Directory structure:
  - `include/` for headers
  - `src/simc_core/` for physics
  - `src/simc_io/` for I/O
  - `data/` for config and matrices
- Library organization:
  - `simc_core` (core physics)
  - `simc_io` (ROOT I/O)
- Basic compilation infrastructure
- Git repository structure

### Files
- `CMakeLists.txt` (top-level)
- `src/simc_core/CMakeLists.txt`
- `src/simc_io/CMakeLists.txt`
- `src/main.cpp`

---

## Version Numbering

Format: `MAJOR.MINOR.PATCH`
- **MAJOR**: Complete phase milestones (e.g., 1.0.0 = all spectrometers)
- **MINOR**: Phase completion (e.g., 0.5.0 = Phase 5)
- **PATCH**: Sub-phase or bug fixes (e.g., 0.5.2 = Phase 5c.2)

Current: **0.5.3** (Phase 5c.3 complete)  
Next: **0.5.4** (Phase 5c.4 - HMS)  
Goal: **1.0.0** (All spectrometers + full validation)

---

## Links
- [Repository](https://github.com/yourusername/simc_cpp)
- [Original Fortran SIMC](https://github.com/ayerbeg/simc_gfortran)
- [Documentation](docs/PHASE_STATUS.md)
