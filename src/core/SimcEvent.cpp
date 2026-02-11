// src/core/SimcEvent.cpp
// Implementation of SIMC event structure

#include "simc/core/SimcEvent.h"
#include <iostream>
#include <iomanip>

ClassImp(simc::SimcEvent)

namespace simc {

// ============================================================================
// Constructor
// ============================================================================
SimcEvent::SimcEvent() {
    Clear();
}

// ============================================================================
// Clear Method
// ============================================================================
void SimcEvent::Clear(Option_t* /*option*/) {
    // Beam and kinematic quantities
    Ein = nu = q = Q2 = W = xbj = epsilon = 0.0;
    
    // Missing quantities
    Em = Pm = Emiss = Pmiss = 0.0;
    Pmx = Pmy = Pmz = 0.0;
    PmPar = PmPer = PmOop = 0.0;
    
    // Recoil system
    Mrec = Trec = 0.0;
    
    // Angles
    theta_pq = phi_pq = theta_tarq = phi_targ = 0.0;
    beta = phi_s = phi_c = 0.0;
    
    // Semi-inclusive
    zhad = pt2 = 0.0;
    
    // Electron arm
    e_delta = e_xptar = e_yptar = e_z = 0.0;
    e_theta = e_phi = e_E = e_P = 0.0;
    
    // Hadron arm
    p_delta = p_xptar = p_yptar = p_z = 0.0;
    p_theta = p_phi = p_E = p_P = 0.0;
    
    // Unit vectors
    ue_x = ue_y = ue_z = 0.0;
    up_x = up_y = up_z = 0.0;
    uq_x = uq_y = uq_z = 0.0;
}

// ============================================================================
// Print Method
// ============================================================================
void SimcEvent::Print(Option_t* /*option*/) const {
    std::cout << "\n========== SIMC Event ==========" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    
    // Beam and kinematics
    std::cout << "Beam Energy:  " << Ein << " MeV" << std::endl;
    std::cout << "Q2:           " << Q2/1e6 << " (GeV/c)^2" << std::endl;
    std::cout << "W:            " << W << " MeV" << std::endl;
    std::cout << "nu:           " << nu << " MeV" << std::endl;
    std::cout << "xbj:          " << xbj << std::endl;
    std::cout << "epsilon:      " << epsilon << std::endl;
    
    // Electron
    std::cout << "\nElectron:" << std::endl;
    std::cout << "  E:      " << e_E << " MeV" << std::endl;
    std::cout << "  P:      " << e_P << " MeV/c" << std::endl;
    std::cout << "  theta:  " << e_theta * constants::degrad << " deg" << std::endl;
    std::cout << "  phi:    " << e_phi * constants::degrad << " deg" << std::endl;
    std::cout << "  delta:  " << e_delta << " %" << std::endl;
    std::cout << "  xptar:  " << e_xptar * 1000.0 << " mrad" << std::endl;
    std::cout << "  yptar:  " << e_yptar * 1000.0 << " mrad" << std::endl;
    
    // Hadron
    std::cout << "\nHadron:" << std::endl;
    std::cout << "  E:      " << p_E << " MeV" << std::endl;
    std::cout << "  P:      " << p_P << " MeV/c" << std::endl;
    std::cout << "  theta:  " << p_theta * constants::degrad << " deg" << std::endl;
    std::cout << "  phi:    " << p_phi * constants::degrad << " deg" << std::endl;
    std::cout << "  delta:  " << p_delta << " %" << std::endl;
    std::cout << "  xptar:  " << p_xptar * 1000.0 << " mrad" << std::endl;
    std::cout << "  yptar:  " << p_yptar * 1000.0 << " mrad" << std::endl;
    
    // Missing quantities
    std::cout << "\nMissing:" << std::endl;
    std::cout << "  Em:     " << Em << " MeV" << std::endl;
    std::cout << "  Pm:     " << Pm << " MeV/c" << std::endl;
    std::cout << "  Pmx:    " << Pmx << " MeV/c" << std::endl;
    std::cout << "  Pmy:    " << Pmy << " MeV/c" << std::endl;
    std::cout << "  Pmz:    " << Pmz << " MeV/c" << std::endl;
    
    std::cout << "================================\n" << std::endl;
}

// ============================================================================
// Set Methods
// ============================================================================
void SimcEvent::SetElectron(const ArmState& arm) {
    e_delta = arm.delta;
    e_xptar = arm.xptar;
    e_yptar = arm.yptar;
    e_z = arm.z;
    e_theta = arm.theta;
    e_phi = arm.phi;
    e_E = arm.E;
    e_P = arm.P;
}

void SimcEvent::SetHadron(const ArmState& arm) {
    p_delta = arm.delta;
    p_xptar = arm.xptar;
    p_yptar = arm.yptar;
    p_z = arm.z;
    p_theta = arm.theta;
    p_phi = arm.phi;
    p_E = arm.E;
    p_P = arm.P;
}

void SimcEvent::SetUnitVectors(const Vector3D& ue, const Vector3D& up, const Vector3D& uq) {
    ue_x = ue.x; ue_y = ue.y; ue_z = ue.z;
    up_x = up.x; up_y = up.y; up_z = up.z;
    uq_x = uq.x; uq_y = uq.y; uq_z = uq.z;
}

// ============================================================================
// Get Methods
// ============================================================================
ArmState SimcEvent::GetElectronState() const {
    ArmState arm;
    arm.delta = e_delta;
    arm.xptar = e_xptar;
    arm.yptar = e_yptar;
    arm.z = e_z;
    arm.theta = e_theta;
    arm.phi = e_phi;
    arm.E = e_E;
    arm.P = e_P;
    return arm;
}

ArmState SimcEvent::GetHadronState() const {
    ArmState arm;
    arm.delta = p_delta;
    arm.xptar = p_xptar;
    arm.yptar = p_yptar;
    arm.z = p_z;
    arm.theta = p_theta;
    arm.phi = p_phi;
    arm.E = p_E;
    arm.P = p_P;
    return arm;
}

Vector3D SimcEvent::GetElectronVector() const {
    return {ue_x, ue_y, ue_z};
}

Vector3D SimcEvent::GetHadronVector() const {
    return {up_x, up_y, up_z};
}

Vector3D SimcEvent::GetPhotonVector() const {
    return {uq_x, uq_y, uq_z};
}

} // namespace simc
