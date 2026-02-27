# Changelog - SIMC C++/ROOT Port

All notable changes to the SIMC C++/ROOT port.

---

## [0.5.4] - 2026-02-13 ✅ HMS Implementation

### Added
- **HMS Spectrometer** (Electron Arm)
  - Complete class implementation (330 lines header, 750 lines source)
  - 14 forward transformations, 3 reconstruction
  - Complex 6-region dipole (Niculescu model with ±6° tilt)
  - 3 quadrupoles with aperture checks
  - HMS-100 collimator (optional octagonal)
  - 3 post-dipole vacuum pipes
  
- **Two-Arm Coincidence**
  - HMS + SHMS working together
  - Proper AND logic validated against Fortran
  - Early exit optimization
  
- **Matrix Parser Improvements**
  - Handles variable coefficients (4 or 5)
  - Fortran-compatible spacing

### Results
- HMS: 83.8% acceptance
- SHMS: 18.7% (of HMS-passed)
- Coincidence: 15.7% ✓

---

## [0.5.3] - 2026-02-11 - Acceptance Generation
- Fixed infinite loop
- Acceptance-based generation
- 100% efficiency, 17.1% SHMS acceptance

---

## [0.5.2] - 2026-02-11 - SHMS Implementation  
- Complete SHMS class (32 apertures)
- 72.5% standalone, 17.1% integrated

---

## Earlier Phases
- 0.5.0-0.5.1: Transport, Cross Sections, Multiple Scattering
- 0.4.0: Radiative Corrections
- 0.3.0: Configuration
- 0.2.0: Data Structures
- 0.1.0: Project Setup

---

## Next: [0.5.5] SOS Implementation (In Progress)
