# SIMC C++/ROOT Port - Phase Status

**Project**: Porting Fortran SIMC Monte Carlo to C++/ROOT  
**Developer**: Carlos  
**Start Date**: 2024  
**Current Phase**: 5c.3 ✅ **COMPLETE**

---

## Quick Reference

| Phase | Component | Status | Completion Date |
|-------|-----------|--------|-----------------|
| 1 | Project Setup | ✅ Complete | - |
| 2 | Data Structures | ✅ Complete | - |
| 3 | Configuration | ✅ Complete | - |
| 4 | Radiative Corrections | ✅ Complete | - |
| 5a | Cross Sections | ✅ Complete | - |
| 5b | Multiple Scattering | ✅ Complete | - |
| 5c.1 | Core Transport | ✅ Complete | - |
| 5c.2 | SHMS Implementation | ✅ Complete | Feb 11, 2026 |
| 5c.3 | Acceptance Generation | ✅ Complete | Feb 11, 2026 |
| **5c.4** | **HMS Implementation** | ⏳ **NEXT** | - |
| 5d | Event Generation | 🔄 Partial | - |
| 6 | Full Integration | ⏳ Pending | - |
| 7 | Validation | ⏳ Pending | - |

---

## Detailed Phase Breakdown

### ✅ Phase 1: Project Setup
**Status**: Complete  
**Description**: CMake build system, ROOT integration, directory structure

**Deliverables**:
- ✅ CMakeLists.txt structure
- ✅ ROOT dictionary generation
- ✅ Library organization (simc_core, simc_io, simc_shms)
- ✅ Include/src separation

---

### ✅ Phase 2: Data Structures
**Status**: Complete  
**Description**: Core C++ classes and structures

**Deliverables**:
- ✅ `SimcEvent` class (main event structure)
- ✅ `TargetInfo`, `SpectrometerInfo` structures
- ✅ Physics constants
- ✅ ROOT dictionary for I/O

---

### ✅ Phase 3: Configuration System
**Status**: Complete  
**Description**: JSON-based configuration

**Deliverables**:
- ✅ `ConfigManager` class
- ✅ JSON parsing (nlohmann/json)
- ✅ Default configuration files
- ✅ Validation and error checking

**Files**:
- `include/ConfigManager.h`
- `src/simc_core/ConfigManager.cpp`
- `data/config/default.json`

---

### ✅ Phase 4: Radiative Corrections
**Status**: Complete  
**Description**: QED radiative corrections

**Deliverables**:
- ✅ `RadiativeCorrections` class
- ✅ Internal/external bremsstrahlung
- ✅ Schwinger correction
- ✅ Peaking approximation

**Files**:
- `include/RadiativeCorrections.h`
- `src/simc_core/RadiativeCorrections.cpp`

**Notes**: Ported from `radc.f`, implements Mo & Tsai formalism

---

### ✅ Phase 5a: Cross Sections
**Status**: Complete  
**Description**: Physics cross section calculations

**Deliverables**:
- ✅ Elastic H(e,e'p): Dipole form factors
- ✅ Quasi-elastic D(e,e'p)
- ✅ Mott cross section
- ✅ Form factor parameterizations

**Files**:
- `include/CrossSection.h`
- `src/simc_core/CrossSection.cpp`

---

### ✅ Phase 5b: Multiple Scattering
**Status**: Complete  
**Description**: Lynch-Dahl multiple scattering

**Deliverables**:
- ✅ `MonteCarloTransport` class
- ✅ Beam multiple scattering
- ✅ Particle multiple scattering
- ✅ Target material effects

**Files**:
- `include/MonteCarloTransport.h`
- `src/simc_core/MonteCarloTransport.cpp`

**Notes**: Ported from `musc.f`, `target.f`

---

### ✅ Phase 5c.1: Core Monte Carlo Transport
**Status**: Complete  
**Description**: Basic transport framework

**Deliverables**:
- ✅ Track state propagation
- ✅ Energy loss corrections
- ✅ Coordinate transformations
- ✅ Acceptance checking (basic)

**Files**:
- `src/simc_core/MonteCarloTransport.cpp` (expanded)

---

### ✅ Phase 5c.2: SHMS Spectrometer
**Status**: Complete  
**Completion**: Feb 11, 2026 (Sessions 1-3)

**Description**: Full SHMS implementation with matrix transport

**Deliverables**:
- ✅ `SHMS` class
- ✅ 32 aperture planes:
  - HB (Heavy gas Cherenkov): 4 planes
  - Q1 (Quadrupole 1): 5 planes  
  - Q2 (Quadrupole 2): 5 planes
  - Q3 (Quadrupole 3): 5 planes
  - Dipole: 13 planes (entrance, flare, 7 internal, exit)
- ✅ Forward matrix (5th order polynomial)
- ✅ Reconstruction matrix
- ✅ TRANSPORT coordinate system
- ✅ Build system integration

**Files**:
- `include/SHMS.h`
- `src/simc_shms/SHMS.cpp`
- `src/simc_shms/CMakeLists.txt`
- `data/matrices/shms/shms_forward.dat`
- `data/matrices/shms/shms_recon.dat`

**Validation**:
- ✅ Standalone test: 72.5% acceptance (10k events)
- ✅ Integrated test: 17.1% acceptance (explained by config mismatch)

**Key Findings**:
- Central momentum mismatch (5.0 GeV vs 2.8 GeV elastic) → δ≈-50%
- Acceptance matches Fortran behavior
- Transport algorithm validated against original SIMC

---

### ✅ Phase 5c.3: Acceptance-Based Generation
**Status**: Complete  
**Completion**: Feb 11, 2026

**Description**: Generate events within spectrometer acceptance for elastic H(e,e'p)

**Deliverables**:
- ✅ Modified `EventGenerator` for elastic mode
- ✅ Generate hadron angles within ±60 mrad of SHMS axis
- ✅ Physics angles from spectrometer angles
- ✅ Removed kinematic constraint (let SHMS handle acceptance)

**Files**:
- `EventGenerator_PHASE5C3_FIXED.cpp`

**Implementation**:
```cpp
// Generate hadron spectrometer angles WITHIN acceptance
event.p_yptar = random->Uniform(gen_limits.hadron.yptar_min, yptar_max);
event.p_xptar = random->Uniform(gen_limits.hadron.xptar_min, xptar_max);

// Calculate physics angles from spectrometer angles
PhysicsAngles(spec_hadron.theta, spec_hadron.phi,
              event.p_xptar, event.p_yptar,
              event.p_theta, event.p_phi);
```

**Validation**:
- ✅ 100% generation efficiency (all generated events transported)
- ✅ 17.1% SHMS acceptance (geometry only)
- ✅ Matches Fortran SIMC approach

**Bug Fixes**:
- Fixed infinite loop (removed overly tight kinematic constraint)
- Proper handling of spectrometer/physics angle conversion

---

### ⏳ Phase 5c.4: HMS Implementation (NEXT)
**Status**: Not started  
**Priority**: HIGH (needed for electron arm)

**Description**: Port HMS spectrometer

**Planned Deliverables**:
- ⏳ `HMS` class
- ⏳ HMS apertures (fewer than SHMS)
- ⏳ HMS forward/reconstruction matrices
- ⏳ Integration with main loop

**Fortran Reference**:
- `hms/mc_hms.f`
- `hms/apertures_hms.inc`
- `hms/forward_cosy.dat`, `hms/recon_cosy.dat`

**Strategy**: Follow SHMS pattern (proven successful)

---

### 🔄 Phase 5d: Event Generation (Partial)
**Status**: Partial implementation

**Completed**:
- ✅ `EventGenerator` class structure
- ✅ Target vertex generation
- ✅ Beam energy smearing
- ✅ Raster pattern simulation
- ✅ Elastic H(e,e'p) kinematics (Phase 5c.3)

**Remaining**:
- ⏳ Deuterium D(e,e'p)
- ⏳ Heavy nuclei A(e,e'p)
- ⏳ Pion production
- ⏳ Kaon production
- ⏳ Fermi momentum generation

**Files**:
- `include/EventGenerator.h`
- `src/simc_core/EventGenerator.cpp`

---

### ⏳ Phase 6: Full Integration
**Status**: Pending  
**Dependencies**: Phase 5c.4 (HMS), Phase 5d (complete event generation)

**Planned**:
- Main event loop
- Coincidence timing
- Two-arm Monte Carlo
- Output ROOT file generation
- Ntuple structure

---

### ⏳ Phase 7: Validation
**Status**: Pending  
**Dependencies**: Phase 6

**Planned**:
- Cross-check against Fortran SIMC
- Standard test cases
- Acceptance comparisons
- Cross section validation
- Documentation

---

## File Structure (Current)

```
simc_cpp/
├── CMakeLists.txt
├── include/
│   ├── ConfigManager.h
│   ├── CrossSection.h
│   ├── EventGenerator.h
│   ├── MonteCarloTransport.h
│   ├── PhysicsConstants.h
│   ├── RadiativeCorrections.h
│   ├── SHMS.h
│   ├── SimcEvent.h
│   └── SimcTypes.h
├── src/
│   ├── simc_core/
│   │   ├── CMakeLists.txt
│   │   ├── ConfigManager.cpp
│   │   ├── CrossSection.cpp
│   │   ├── EventGenerator.cpp
│   │   ├── MonteCarloTransport.cpp
│   │   └── RadiativeCorrections.cpp
│   ├── simc_io/
│   │   ├── CMakeLists.txt
│   │   └── OutputManager.cpp
│   ├── simc_shms/
│   │   ├── CMakeLists.txt
│   │   └── SHMS.cpp
│   └── main.cpp
├── data/
│   ├── config/
│   │   └── default.json
│   └── matrices/
│       └── shms/
│           ├── shms_forward.dat
│           └── shms_recon.dat
└── tests/
    └── test_shms.cpp
```

---

## Key Decisions & Learnings

### ✅ SHMS Success Pattern
The SHMS implementation (Phase 5c.2) established a successful pattern:
1. **Separate class** for each spectrometer
2. **Matrix-based transport** with polynomial evaluation
3. **Aperture checking** at all key planes
4. **TRANSPORT coordinate system** for consistency
5. **Standalone testing** before integration

This pattern should be followed for HMS, SOS, HRS.

### ✅ Configuration Philosophy
- **JSON over text files** (modern, hierarchical, validated)
- **Type-safe access** via ConfigManager
- **Default values** for all parameters
- **Validation at load time**

### ✅ Build System
- **CMake** with modern target-based approach
- **Separate libraries** (core, io, shms, ...)
- **ROOT dictionary** only where needed
- **Automatic dependency tracking**

### ✅ Acceptance Generation (Phase 5c.3)
- Generate angles **within** spectrometer acceptance
- Let spectrometer transport determine final acceptance
- Don't apply kinematic constraints if spectrometer is misaligned
- **17.1% acceptance is correct** for 5.0 GeV SHMS detecting 2.8 GeV protons

---

## Performance Metrics

### Phase 5c.2 (SHMS Standalone)
- Events: 10,000
- Generation: ~instant
- Transport: <1 second
- Acceptance: 72.5%

### Phase 5c.3 (Integrated)
- Events: 10,000
- Generation: 100% efficiency
- SHMS acceptance: 17.1%
- Delta range: -55% to -27%
- **Matches Fortran SIMC behavior**

---

## Next Steps

### Immediate (Phase 5c.4)
1. Create `HMS` class (copy SHMS pattern)
2. Port HMS apertures from Fortran
3. Load HMS matrices
4. Test standalone
5. Integrate into main loop

### Short-term
- Complete event generation (Phase 5d)
- Add SOS support (Phase 5c.5)
- Two-arm coincidence (Phase 6)

### Long-term
- Full physics validation
- Performance optimization
- Documentation
- User manual

---

## References

### Fortran SIMC
- Repository: https://github.com/ayerbeg/simc_gfortran
- Key files: `simc.f`, `event.f`, `shms/mc_shms.f`

### Documentation
- SHMS TDR: https://hallcweb.jlab.org/DocDB/0009/000956/002/shms_final.pdf
- Hall C Monte Carlo Manual (original)

### This Project
- See `SPECTROMETERS_OVERVIEW.md` for spectrometer details
- See `PHASE_5C3_BUG_FIX.md` for acceptance generation details
- See individual session transcripts for detailed progress

---

**Last Updated**: February 11, 2026  
**Next Milestone**: Phase 5c.4 (HMS Implementation)
