# Hall C Spectrometers in SIMC

## Overview

Based on the Fortran SIMC codebase, there are **6 spectrometer systems** supported:

| ID | Name | Description | Status in C++ Port |
|----|------|-------------|-------------------|
| 1  | **HMS** | High Momentum Spectrometer | ⏳ Not yet ported |
| 2  | **SOS** | Short Orbit Spectrometer | ⏳ Not yet ported |
| 3  | **HRSr** | Hall A HRS Right arm | ⏳ Not yet ported |
| 4  | **HRSl** | Hall A HRS Left arm | ⏳ Not yet ported |
| 5/6 | **SHMS** | Super High Momentum Spectrometer | ✅ **COMPLETED Phase 5c.2** |
| 7/8 | **CALO** | Calorimeter (photon detection) | ⏳ Not yet ported |

## Spectrometer Details

### 1. HMS (High Momentum Spectrometer)
- **Hall**: Hall C
- **Fortran files**: `hms/mc_hms.f`, `hms/mc_hms_recon.f`, `hms/mc_hms_hut.f`
- **Matrix files**: `hms/forward_cosy.dat`, `hms/recon_cosy.dat`
- **Apertures**: `hms/apertures_hms.inc`
- **Features**:
  - 3 quadrupoles + 1 dipole
  - Typical central momentum: 0.4-7.5 GeV/c
  - Solid angle: ~6 msr
  - Momentum acceptance: ±10%
  - Angular acceptance: ±28 mrad (horizontal), ±50 mrad (vertical)

### 2. SOS (Short Orbit Spectrometer)
- **Hall**: Hall C
- **Fortran files**: `sos/mc_sos.f`, `sos/mc_sos_recon.f`, `sos/mc_sos_hut.f`
- **Matrix files**: `sos/forward_cosy.dat`, `sos/recon_cosy.dat`
- **Apertures**: `sos/apertures_sos.inc`
- **Features**:
  - 1 quadrupole + 2 dipoles
  - Typical central momentum: 0.17-1.7 GeV/c
  - Solid angle: ~7 msr
  - Momentum acceptance: ±40%
  - Large acceptance spectrometer

### 3. HRSr (Hall A HRS Right)
- **Hall**: Hall A
- **Fortran files**: `hrsr/mc_hrsr.f`, `hrsr/mc_hrsr_recon.f`, `hrsr/mc_hrsr_hut.f`
- **Matrix files**: `hrsr/hrs_forward_cosy.dat`, `hrsr/hrs_recon_cosy.dat`
- **Apertures**: `hrsr/apertures_hrsr.inc`
- **Features**:
  - 3 quadrupoles + 1 dipole
  - High resolution spectrometer
  - Momentum range: 0.3-4.0 GeV/c
  - Δp/p ~ 10^-4

### 4. HRSl (Hall A HRS Left)
- **Hall**: Hall A
- **Fortran files**: `hrsl/mc_hrsl.f`, `hrsl/mc_hrsl_recon.f`, `hrsl/mc_hrsl_hut.f`
- **Matrix files**: `hrsl/hrs_forward_cosy.dat`, `hrsl/hrs_recon_cosy.dat`
- **Apertures**: `hrsl/apertures_hrsl.inc`
- **Features**:
  - Mirror image of HRSr
  - Same specifications as HRSr

### 5/6. SHMS (Super High Momentum Spectrometer)
- **Hall**: Hall C (12 GeV upgrade)
- **Fortran files**: `shms/mc_shms.f`, `shms/mc_shms_recon.f`, `shms/mc_shms_hut.f`
- **Matrix files**: `shms/shms_forward.dat`, `shms/shms_recon.dat`
- **Apertures**: `shms/apertures_shms.inc`
- **C++ Implementation**: ✅ **COMPLETE** (Phase 5c.2)
  - Class: `SHMS` in `include/SHMS.h`, `src/simc_shms/SHMS.cpp`
  - Features implemented:
    - 32 aperture checks (HB, Q1, Q2, Q3, Dipole with 7 internal planes)
    - Forward matrix (5th order polynomial)
    - Reconstruction matrix
    - Coordinate rotations (TRANSPORT system)
- **Features**:
  - 3 quadrupoles + 1 dipole + optional heavy gas Cherenkov
  - High momentum: 2-11 GeV/c (can go up to 16 GeV/c)
  - Large acceptance: ±40 mrad (horizontal), ±100 mrad (vertical)
  - Momentum acceptance: -10% to +22%
  - Used for hadron detection at 12 GeV

### 7/8. CALO (Calorimeter)
- **Hall**: Hall C
- **Fortran files**: `calo/mc_calo.f`, `calo/mc_calo_recon.f`
- **Purpose**: Photon detection (π⁰ → γγ decay)
- **Features**:
  - Specialized for neutral particle detection
  - Works with `doing_pizero` flag
  - Tracks both photons from pion decay

## Fortran Code Structure

### Main Monte Carlo Files (per spectrometer)
Each spectrometer has 3 main routines:

1. **mc_XXX.f**: Transport through spectrometer
   - Aperture checking
   - Multiple scattering
   - Matrix application
   
2. **mc_XXX_recon.f**: Reconstruction
   - Focal plane → target coordinates
   - Resolution smearing
   
3. **mc_XXX_hut.f**: Detector hut
   - Detector acceptances
   - Wire chamber hits
   - Particle identification

### Configuration (from simc.f)
```fortran
! Spectrometer selection (electron_arm, hadron_arm):
! 1 = HMS
! 2 = SOS  
! 3 = HRSr
! 4 = HRSl
! 5 = SHMS (large angle mode)
! 6 = SHMS (small angle mode)
! 7 = CALO (right)
! 8 = CALO (left)
```

### Monte Carlo Call Pattern (from simc.f, lines ~300-360)
```fortran
if (hadron_arm.eq.1) then
  call mc_hms(...)
else if (hadron_arm.eq.2) then
  call mc_sos(...)
else if (hadron_arm.eq.3) then
  call mc_hrsr(...)
else if (hadron_arm.eq.4) then
  call mc_hrsl(...)
else if (hadron_arm.eq.5 .or. hadron_arm.eq.6) then
  call mc_shms(..., hadron_arm, ...)
else if (hadron_arm.eq.7 .or. hadron_arm.eq.8) then
  call mc_calo(..., hadron_arm, ...)
endif
```

## C++ Port Status

### ✅ Completed: SHMS
- **Phase 5c.2** (Session 1-3): Full SHMS implementation
  - 32 aperture planes
  - 5th order forward matrix
  - Reconstruction matrix
  - Validated: 72.5% standalone acceptance, 17.1% integrated

### ⏳ Next Priorities

**Recommended order:**

1. **HMS** (Phase 5c.4)
   - Most commonly used for electron arm
   - Simpler than SHMS (fewer apertures)
   - Would enable full (e,e'p) coincidence simulation
   
2. **SOS** (Phase 5c.5)
   - Complementary to HMS
   - Large acceptance
   
3. **HRS (l/r)** (Phase 5c.6-7)
   - Hall A support
   - High resolution capabilities
   
4. **CALO** (Phase 5c.8)
   - Specialized, lower priority
   - Only needed for π⁰ studies

## Implementation Template

Based on SHMS success, each spectrometer port should follow this pattern:

```cpp
class HMS {
private:
    struct Aperture {
        double x_min, x_max;
        double y_min, y_max;
    };
    std::vector<Aperture> apertures_;
    
    // Matrix elements
    std::vector<TransformationClass> forward_matrix_;
    std::vector<TransformationClass> recon_matrix_;
    
public:
    bool Transport(const TrackState& initial, 
                   TrackState& final,
                   bool& accepted);
    
    void CalculateReconstructed(const TrackState& focal_plane,
                                TrackState& target);
};
```

## Matrix Files

Each spectrometer needs:
- **Forward matrix**: Target → Focal Plane
- **Reconstruction matrix**: Focal Plane → Target

Format is COSY INFINITY output (polynomial coefficients).

## References

- SHMS TDR: https://hallcweb.jlab.org/DocDB/0009/000956/002/shms_final.pdf
- HMS/SOS: Hall C website documentation
- HRS: Hall A documentation

## Next Steps for C++ Port

1. Create `HMS` class following SHMS pattern
2. Port HMS apertures from `hms/apertures_hms.inc`
3. Load HMS matrix files
4. Integrate into main simulation loop
5. Test with electron arm configuration
