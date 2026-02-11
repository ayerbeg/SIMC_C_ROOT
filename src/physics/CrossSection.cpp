// src/physics/CrossSection.cpp
// COMPLETE PORT from physics_proton.f, physics_pion.f, physics_kaon.f
// 
// PRIMARY SOURCES:
// - physics_proton.f: sigep(), deForest(), fofa_best_fit(), sigMott()
// - physics_pion.f: pion electroproduction (TODO)
// - physics_kaon.f: kaon electroproduction (TODO)

#include "simc/physics/CrossSection.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/SimcConstants.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace simc {
using namespace constants;

// ============================================================================
// ElasticCrossSection Implementation
// PORT FROM: physics_proton.f
// ============================================================================

double ElasticCrossSection::Calculate(const SimcEvent& evt) const {
    // PORT FROM: real*8 function sigep(vertex) in physics_proton.f
    // elastic cross section, units are set by sigMott.f (microbarn/sr)
    
    double q4sq = evt.Q2;  // Q2 is already positive in our convention
    
    // Get form factors - note Fortran passes -q4sq/hbarc**2
    double GE, GM;
    FofaBestFit(-q4sq / (hbarc * hbarc), GE, GM);
    
    // Calculate structure functions
    double qmu4mp = q4sq / (4.0 * Mp2);
    double W1p = GM * GM * qmu4mp;
    double W2p = (GE*GE + GM*GM*qmu4mp) / (1.0 + qmu4mp);
    
    // Rosenbluth structure function
    double tan_half = std::tan(evt.e_theta / 2.0);
    double Wp = W2p + 2.0 * W1p * tan_half * tan_half;
    
    // Mott cross section
    double sig_mott = SigMott(evt.e_E, evt.e_theta, evt.Q2);
    
    // Recoil factor
    double recoil = evt.e_E / evt.Ein;
    
    // Total elastic cross section
    return sig_mott * recoil * Wp;  // microbarn/sr
}

bool ElasticCrossSection::IsPhysical(const SimcEvent& evt) const {
    // Check if elastic kinematics are physical
    
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0) return false;
    if (evt.e_E >= evt.Ein) return false;
    if (evt.e_theta <= 0.0 || evt.e_theta >= pi) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // For elastic: W should equal Mp (allow 100 MeV tolerance)
    if (std::abs(evt.W - Mp) > 100.0) return false;
    
    return true;
}

double ElasticCrossSection::SigMott(double e0, double theta, double Q2) {
    // PORT FROM: real*8 function sigMott(e0,theta,Q2) in physics_proton.f
    // The Mott cross section (for a point nucleus) in microbarns/sr.
    
    double sig = (2.0 * alpha * hbarc * e0 * std::cos(theta / 2.0) / Q2);
    sig = sig * sig;
    
    return sig * 1.0e4;  // fm**2 --> microbarns
}

void ElasticCrossSection::FofaBestFit(double qsquar, double& GE, double& GM) {
    // PORT FROM: subroutine fofa_best_fit(qsquar,GE,GM) in physics_proton.f
    //
    // csa 9/14/98 -- This calculates the form factors Gep and Gmp using
    // Peter Bosted's fit to world data (Phys. Rev. C 51, 409, Eqs. 4
    // and 5 or, alternatively, Eqs. 6)
    //
    // INPUT:  qsquar = -Q^2/(hbarc^2) in fm^-2
    // OUTPUT: GE, GM = electric and magnetic form factors (dimensionless)
    
    constexpr double mu_p = 2.793;  // Proton magnetic moment
    
    // Q2 in GeV^2
    double Q2 = -qsquar * (hbarc * hbarc) * 1.0e-6;
    double Q = std::sqrt(std::max(Q2, 0.0));
    
    double Q2_pow = Q * Q;
    double Q3 = Q2_pow * Q;
    double Q4 = Q2_pow * Q2_pow;
    double Q5 = Q4 * Q;
    
    // Use Eqs 4, 5 from Bosted paper:
    double denom = 1.0 + 0.62*Q + 0.68*Q2_pow + 2.8*Q3 + 0.83*Q4;
    GE = 1.0 / denom;
    
    denom = 1.0 + 0.35*Q + 2.44*Q2_pow + 0.5*Q3 + 1.04*Q4 + 0.34*Q5;
    GM = mu_p / denom;
    
    // Alternative Eqs 6 (commented out in Fortran):
    // denom = 1.0 + 0.14*Q + 3.01*Q2 + 0.02*Q3 + 1.20*Q4 + 0.32*Q5;
    // GE = 1.0 / denom;
    // GM = mu_p / denom;
}

// ============================================================================
// QuasiElasticCrossSection Implementation
// PORT FROM: physics_proton.f
// ============================================================================

QuasiElasticCrossSection::QuasiElasticCrossSection(int deforest_flag)
    : deforest_flag_(deforest_flag) {
}

double QuasiElasticCrossSection::Calculate(const SimcEvent& evt) const {
    // PORT FROM: real*8 function deForest(ev) in physics_proton.f
    //
    // Compute deForest sigcc cross-section, according to value of DEFOREST_FLAG:
    //   Flag = 0    -- use sigcc1
    //   Flag = 1    -- use sigcc2
    //   Flag = -1   -- use sigcc1 ONSHELL, replacing Ebar with E = E'-nu
    //                  and qbar with q (4-vector)
    //
    // N.B. Beware of deForest's metric when taking all those 4-vector inner
    // products in sigcc2 ... it is (-1,1,1,1)!
    // Here, I've defined all the inner products with the regular signs,
    // and then put them in the structure function
    // formulas with reversed signs compared to deForest.
    //
    // UNITS: microbarn*MeV^2/sr^2
    // (need to multiply by spectral function S(E,p) in MeV^-4 to get full 6-fold cross section)
    
    double q4sq = -evt.Q2;  // Fortran uses negative Q2
    double q2 = evt.q * evt.q;
    
    double ebar, qbsq;
    if (deforest_flag_ >= 0) {
        ebar = std::sqrt(evt.Pm*evt.Pm + Mp2);  // Mh2 = Mp2 for protons
        qbsq = (evt.p_E - ebar)*(evt.p_E - ebar) - q2;
    } else {
        ebar = evt.p_E - evt.nu;
        qbsq = q4sq;
    }
    
    // Calculate sin(gamma) where gamma is angle between q and p
    // sin_gamma = sin(angle between q and p)
    double cos_gamma = evt.uq_x*evt.up_x + evt.uq_y*evt.up_y + evt.uq_z*evt.up_z;
    double sin_gamma = 1.0 - cos_gamma * cos_gamma;
    
    if (sin_gamma < 0.0) {
        std::cerr << "WARNING: deForest came up with sin_gamma = " << sin_gamma << std::endl;
        sin_gamma = 0.0;
    }
    sin_gamma = std::sqrt(sin_gamma);
    
    // Calculate cos(phi) where phi is out-of-plane angle
    double cos_phi = 0.0;
    if (sin_gamma != 0.0) {
        cos_phi = (evt.uq_y*(evt.uq_y*evt.up_z - evt.uq_z*evt.up_y)
                  - evt.uq_x*(evt.uq_z*evt.up_x - evt.uq_x*evt.up_z))
                  / sin_gamma / std::sqrt(1.0 - evt.uq_z*evt.uq_z);
    }
    
    if (std::abs(cos_phi) > 1.0) {
        // Set to +/-1, warn if > 1.e-10
        if ((std::abs(cos_phi) - 1.0) > 1.0e-10) {
            std::cerr << "WARNING: deForest give cos_phi = " << cos_phi << std::endl;
        }
        cos_phi = (cos_phi > 0.0) ? 1.0 : -1.0;
    }
    
    // Get form factors
    double GE, GM;
    ElasticCrossSection::FofaBestFit(q4sq / (hbarc * hbarc), GE, GM);
    
    // Calculate Pauli and Dirac form factors
    double qmu4mp = q4sq / (4.0 * Mp2);
    double f1 = (GE - GM*qmu4mp) / (1.0 - qmu4mp);
    double kf2 = (GM - GE) / (1.0 - qmu4mp);
    double f1sq = f1 * f1;
    double kf2_over_2m_allsq = kf2*kf2 / (4.0 * Mp2);
    
    // Calculate kinematic terms
    double tan_half = std::tan(evt.e_theta / 2.0);
    double termC = (q4sq / q2) * (q4sq / q2);
    double termT = tan_half*tan_half - q4sq / (2.0 * q2);
    double termS = tan_half*tan_half - (q4sq / q2) * cos_phi*cos_phi;
    double termI = (-q4sq / q2) * std::sqrt(tan_half*tan_half - q4sq/q2) * cos_phi;
    
    // Calculate structure functions W_C, W_T, W_S, W_I
    double WC, WT, WS, WI;
    
    if (deforest_flag_ <= 0) {
        // Use sigcc1 (or sigcc1 ONSHELL)
        double sumFF1 = (f1 + kf2) * (f1 + kf2);
        double sumFF2 = f1sq - qbsq*kf2*kf2 / (4.0 * Mp2);
        
        WC = ((ebar + evt.p_E)*(ebar + evt.p_E)) * sumFF2 - q2 * sumFF1;
        WT = -2.0 * qbsq * sumFF1;
        WS = 4.0 * (evt.p_P*evt.p_P) * (sin_gamma*sin_gamma) * sumFF2;
        WI = -4.0 * (ebar + evt.p_E) * evt.p_P * sin_gamma * sumFF2;
    } else {
        // Use sigcc2
        double pbarp = ebar*evt.p_E - evt.p_P*(evt.up_x*evt.Pmx + evt.up_y*evt.Pmy + evt.up_z*evt.Pmz);
        double pbarq = ebar*evt.nu - evt.q*(evt.uq_x*evt.Pmx + evt.uq_y*evt.Pmy + evt.uq_z*evt.Pmz);
        double pq = evt.p_E*evt.nu - evt.p_P * evt.q * (evt.up_x*evt.uq_x + evt.up_y*evt.uq_y + evt.up_z*evt.uq_z);
        double qbarq = (evt.p_E - ebar)*evt.nu - q2;
        
        WC = (ebar*evt.p_E + (-pbarp + Mp2)/2.0) * f1sq - q2*f1*kf2/2.0
             - ((-pbarq*evt.p_E - pq*ebar)*evt.nu + ebar*evt.p_E*q4sq
                + pbarq*pq - (-pbarp - Mp2)/2.0*q2) * kf2_over_2m_allsq;
        
        WT = -(-pbarp + Mp2)*f1sq - qbarq*f1*kf2
             + (2.0*pbarq*pq + (-pbarp - Mp2)*q4sq) * kf2_over_2m_allsq;
        
        WS = (evt.p_P*sin_gamma)*(evt.p_P*sin_gamma) * (f1sq - q4sq*kf2_over_2m_allsq);
        
        WI = evt.p_P*sin_gamma * (-(ebar + evt.p_E)*f1sq
                                  + ((-pbarq - pq)*evt.nu + (ebar + evt.p_E)*q4sq)*kf2_over_2m_allsq);
    }
    
    // Sum all structure function contributions
    double allsum = termC*WC + termT*WT + termS*WS + termI*WI;
    
    if (deforest_flag_ <= 0) {
        allsum = allsum / 4.0;
    }
    
    // Calculate final cross section
    double deforest_xs = ElasticCrossSection::SigMott(evt.e_E, evt.e_theta, evt.Q2)
                         * evt.p_P * allsum / ebar;
    
    return deforest_xs;  // microbarn*MeV^2/sr^2
}

bool QuasiElasticCrossSection::IsPhysical(const SimcEvent& evt) const {
    // Check quasi-elastic kinematics
    
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // Missing momentum should be reasonable (< 1 GeV/c)
    if (evt.Pm < 0.0 || evt.Pm > 1000.0) return false;
    
    // Missing energy range
    // For light nuclei: Em ~ -50 to +50 MeV (binding + recoil)
    // For high-Q2 scattering: Em can be several hundred MeV (recoil system KE)
    // Allow wide range: -500 to +1000 MeV
    if (evt.Em < -500.0 || evt.Em > 1000.0) return false;
    
    return true;
}

// ============================================================================
// PionCrossSection Implementation
// PORT FROM: physics_pion.f (TODO - needs MAID tables)
// ============================================================================

PionCrossSection::PionCrossSection(PionType type) 
    : pion_type_(type) {
}

double PionCrossSection::Calculate(const SimcEvent& evt) const {
    // PORT FROM: physics_pion.f
    // TODO: Full implementation requires MAID table integration
    // This is a placeholder - full implementation would read MAID tables
    
    double Q2_GeV2 = evt.Q2 / 1.0e6;
    double W_GeV = evt.W / 1.0e3;
    
    // Simple delta resonance peak
    double W_delta = 1.232;  // GeV
    double Gamma = 0.117;    // GeV
    
    // Breit-Wigner resonance
    double denominator = (W_GeV - W_delta)*(W_GeV - W_delta) + Gamma*Gamma/4.0;
    double resonance = Gamma*Gamma / (4.0 * denominator);
    
    // Q2 dependence (simple dipole)
    double Lambda2 = 0.71 * 0.71;  // GeV^2
    double FF = 1.0 / ((1.0 + Q2_GeV2/Lambda2) * (1.0 + Q2_GeV2/Lambda2));
    
    // Rough normalization (microbarns/sr^2/MeV)
    double norm = 100.0;
    
    return norm * resonance * FF;
}

bool PionCrossSection::IsPhysical(const SimcEvent& evt) const {
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // W must be above pion production threshold
    double W_threshold = Mp + Mpi;  // ~ 1078 MeV
    if (evt.W < W_threshold) return false;
    if (evt.W > 3000.0) return false;  // 3 GeV
    
    return true;
}

double PionCrossSection::HadronicTensor(const SimcEvent& /* evt */) const {
    // TODO: Implement hadronic tensor for pion production
    return 1.0;
}

// ============================================================================
// KaonCrossSection Implementation
// PORT FROM: physics_kaon.f (TODO - needs Saghai model)
// ============================================================================

KaonCrossSection::KaonCrossSection(KaonType type)
    : kaon_type_(type) {
}

double KaonCrossSection::Calculate(const SimcEvent& evt) const {
    // PORT FROM: physics_kaon.f
    // TODO: Full implementation requires Saghai model
    
    double W_GeV = evt.W / 1.0e3;
    
    // Kaon production threshold
    double W_threshold = (Mp + Mk) / 1.0e3;  // GeV
    if (W_GeV < W_threshold) return 0.0;
    
    // Simple Regge model parameterization
    double t = evt.Q2 / 1.0e6;  // -t in GeV^2
    double s = W_GeV * W_GeV;   // s in GeV^2
    
    // Regge trajectory
    double alpha_prime = 0.9;  // GeV^-2
    double alpha_0 = 0.5;
    double alpha_t = alpha_0 + alpha_prime * t;
    
    // Cross section ~ s^(2*alpha-2)
    double xs = std::pow(s, 2.0*alpha_t - 2.0);
    
    // Normalization (microbarns)
    double norm = 10.0;
    
    return norm * xs;
}

bool KaonCrossSection::IsPhysical(const SimcEvent& evt) const {
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // W must be above kaon production threshold
    double W_threshold = Mp + Mk;  // ~ 1432 MeV
    if (evt.W < W_threshold) return false;
    if (evt.W > 5000.0) return false;  // 5 GeV
    
    return true;
}

// ============================================================================
// CrossSectionFactory Implementation
// ============================================================================

std::unique_ptr<CrossSectionBase> CrossSectionFactory::Create(ReactionType type) {
    switch (type) {
        case ReactionType::ELASTIC:
            return std::make_unique<ElasticCrossSection>();
            
        case ReactionType::QUASIELASTIC:
            return std::make_unique<QuasiElasticCrossSection>(0);  // Default flag = 0
            
        case ReactionType::PION_PRODUCTION:
            return std::make_unique<PionCrossSection>(PionType::PI_PLUS);
            
        case ReactionType::KAON_PRODUCTION:
            return std::make_unique<KaonCrossSection>(KaonType::K_PLUS);
            
        default:
            throw std::runtime_error("Unsupported reaction type");
    }
}

std::unique_ptr<CrossSectionBase> CrossSectionFactory::CreateFromConfig(
    const ConfigManager& config) {
    
    ReactionType type = config.GetReactionType();
    auto xs = Create(type);
    
    // Configure specific parameters based on reaction
    // TODO: Add pion type, kaon type, deForest flag from config
    
    return xs;
}

// ============================================================================
// RadiativeCorrections Implementation
// PORT FROM: radc.f (TODO - full implementation)
// ============================================================================

RadiativeCorrections::RadiativeCorrections(double bt) : bt_(bt) {
}

double RadiativeCorrections::GetCorrectionFactor(double Ein, double Ee, 
                                                 double theta, int Z) const {
    // PORT FROM: radc.f
    // TODO: Full implementation
    
    double internal = InternalCorrection(Ein, Ee, theta);
    double external = ExternalCorrection(Ein, Ee, Z, 0.1);
    
    return (1.0 + internal) * (1.0 + external);
}

double RadiativeCorrections::InternalCorrection(double Ein, double Ee, double theta) const {
    // PORT FROM: radc.f - basicrad()
    // TODO: Full implementation
    
    double Q2 = 2.0 * Ein * Ee * (1.0 - std::cos(theta));
    
    double delta_schwinger = (alpha / pi) * (std::log(Q2 / Me2) - 1.0);
    double delta_vertex = -(alpha / pi) * (pi*pi / 6.0);
    double delta_vacuum = (alpha / (3.0 * pi)) * std::log(Q2 / Me2);
    
    return delta_schwinger + delta_vertex + delta_vacuum;
}

double RadiativeCorrections::ExternalCorrection(double Ein, double Ee, 
                                                int Z, double thickness) const {
    // PORT FROM: brem.f
    // TODO: Full implementation
    
    double X0 = 36.1 / Z;
    double t = thickness / X0;
    double E_avg = (Ein + Ee) / 2.0;
    double delta_E = bt_ * E_avg * t;
    
    return delta_E / Ee;
}

double RadiativeCorrections::RadiativeTail(double Ein, double Ee, double theta) const {
    // TODO: Full implementation
    
    double Q2 = 2.0 * Ein * Ee * (1.0 - std::cos(theta));
    double nu = Ein - Ee;
    double b = (alpha / pi) * (std::log(Q2 / Me2) - 1.0);
    double tail = b * nu / Ein;
    
    return tail;
}

} // namespace simc
