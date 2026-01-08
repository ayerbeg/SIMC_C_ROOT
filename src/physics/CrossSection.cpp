// src/physics/CrossSection.cpp
// Implementation of cross section calculations
// Ported from physics_proton.f, physics_pion.f, physics_kaon.f

#include "simc/CrossSection.h"
#include "simc/ConfigManager.h"
#include "simc/SimcConstants.h"
#include <cmath>
#include <stdexcept>

namespace simc {
using namespace constants;

// ============================================================================
// ElasticCrossSection Implementation
// ============================================================================

double ElasticCrossSection::Calculate(const SimcEvent& evt) const {
    // Port from sigep() in physics_proton.f
    // elastic cross section, units are microbarn/sr
    
    double Q2_GeV2 = evt.Q2 / 1.0e6;  // Convert MeV^2 to GeV^2
    
    // Get form factors
    double GE = GetGE(Q2_GeV2);
    double GM = GetGM(Q2_GeV2);
    
    // Calculate structure functions
    double qmu4mp = evt.Q2 / (4.0 * Mp2);
    double W1p = GM * GM * qmu4mp;
    double W2p = (GE*GE + GM*GM*qmu4mp) / (1.0 + qmu4mp);
    
    // Rosenbluth structure function
    double tan_half = std::tan(evt.e_theta / 2.0);
    double Wp = W2p + 2.0 * W1p * tan_half * tan_half;
    
    // Mott cross section
    double sig_mott = MottCrossSection(evt.Ein, evt.e_theta);
    
    // Recoil factor
    double recoil = evt.e_E / evt.Ein;
    
    // Total elastic cross section
    return sig_mott * recoil * Wp;  // microbarn/sr
}

bool ElasticCrossSection::IsPhysical(const SimcEvent& evt) const {
    // Check if elastic kinematics are physical
    
    // Energy must be positive
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0) return false;
    
    // Must have scattered electron
    if (evt.e_E >= evt.Ein) return false;
    
    // Angle must be reasonable
    if (evt.e_theta <= 0.0 || evt.e_theta >= pi) return false;
    
    // Q2 must be positive
    if (evt.Q2 <= 0.0) return false;
    
    // Check elastic peak constraint
    // For elastic: W should equal Mp
    if (std::abs(evt.W - Mp) > 100.0) return false;  // 100 MeV tolerance
    
    return true;
}

double ElasticCrossSection::MottCrossSection(double Ein, double theta) {
    // Port from sigMott() in physics_proton.f
    // The Mott cross section (for a point nucleus) in microbarns/sr
    
    double Q2 = 2.0 * Ein * Ein * (1.0 - std::cos(theta));
    
    // Mott formula: (2 * alpha * hbarc * E * cos(theta/2) / Q2)^2
    double cos_half = std::cos(theta / 2.0);
    double numerator = 2.0 * alpha * hbarc * Ein * cos_half;
    double sig = (numerator / Q2) * (numerator / Q2);
    
    // Convert fm^2 to microbarns (1 fm^2 = 10 microbarns)
    return sig * 1.0e4;  // microbarn/sr
}

double ElasticCrossSection::RecoilFactor(double Ein, double theta, double M_target) {
    // Recoil factor for elastic scattering
    // Scattered energy: E' = E / (1 + 2E*sin^2(theta/2)/M)
    
    double sin_half = std::sin(theta / 2.0);
    double denominator = 1.0 + (2.0 * Ein * sin_half * sin_half) / M_target;
    
    return 1.0 / denominator;
}

double ElasticCrossSection::GetGE(double Q2) {
    // Port from fofa_best_fit() in physics_proton.f
    // Peter Bosted's fit to world data (Phys. Rev. C 51, 409)
    // Q2 in GeV^2
    
    double Q = std::sqrt(std::max(Q2, 0.0));
    double Q2_pow = Q * Q;
    double Q3 = Q2_pow * Q;
    double Q4 = Q2_pow * Q2_pow;
    
    // Equation 4 from Bosted
    double denom = 1.0 + 0.62*Q + 0.68*Q2_pow + 2.8*Q3 + 0.83*Q4;
    
    return 1.0 / denom;
}

double ElasticCrossSection::GetGM(double Q2) {
    // Port from fofa_best_fit() in physics_proton.f
    // Magnetic form factor
    // Q2 in GeV^2
    
    const double mu_p = 2.793;  // Proton magnetic moment
    
    double Q = std::sqrt(std::max(Q2, 0.0));
    double Q2_pow = Q * Q;
    double Q3 = Q2_pow * Q;
    double Q4 = Q2_pow * Q2_pow;
    double Q5 = Q4 * Q;
    
    // Equation 5 from Bosted
    double denom = 1.0 + 0.35*Q + 2.44*Q2_pow + 0.5*Q3 + 1.04*Q4 + 0.34*Q5;
    
    return mu_p / denom;
}

double ElasticCrossSection::DipoleFF(double Q2, double Lambda) {
    // Simple dipole form factor
    // Q2 in GeV^2, Lambda in GeV
    
    double Lambda2 = Lambda * Lambda;
    double denominator = (1.0 + Q2 / Lambda2);
    
    return 1.0 / (denominator * denominator);
}

// ============================================================================
// QuasiElasticCrossSection Implementation
// ============================================================================

double QuasiElasticCrossSection::Calculate(const SimcEvent& evt) const {
    // Port from deForest() in physics_proton.f
    // Quasi-elastic (e,e'p) cross section with off-shell corrections
    
    // Get elastic cross section as baseline
    ElasticCrossSection elastic;
    double sig_elastic = elastic.Calculate(evt);
    
    // Off-shell correction factor
    double offshell = OffShellCorrection(evt.Pm, evt.Em);
    
    // Apply correction
    return sig_elastic * offshell;
}

bool QuasiElasticCrossSection::IsPhysical(const SimcEvent& evt) const {
    // Check quasi-elastic kinematics
    
    // Basic kinematic checks
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // Missing momentum should be reasonable (< 1 GeV/c)
    if (evt.Pm < 0.0 || evt.Pm > 1000.0) return false;
    
    // Missing energy should be near binding energy
    // For nuclei, typical binding ~ -10 to -50 MeV
    if (evt.Em < -500.0 || evt.Em > 200.0) return false;
    
    return true;
}

double QuasiElasticCrossSection::OffShellCorrection(double Pm, double Em) {
    // Off-shell correction for bound nucleons
    // This is a simplified version - full version would use spectral functions
    
    // For now, use a simple Gaussian fall-off with missing momentum
    // Real implementation would integrate over spectral function
    
    double Pm_fermi = 250.0;  // Fermi momentum ~ 250 MeV/c
    double correction = std::exp(-Pm*Pm / (2.0 * Pm_fermi * Pm_fermi));
    
    // Binding energy correction
    double E_binding = -8.0;  // Typical nuclear binding ~ 8 MeV
    double E_width = 20.0;
    correction *= std::exp(-(Em - E_binding)*(Em - E_binding) / (2.0 * E_width * E_width));
    
    return correction;
}

// ============================================================================
// PionCrossSection Implementation
// ============================================================================

PionCrossSection::PionCrossSection(PionType type) 
    : pion_type_(type) {
}

double PionCrossSection::Calculate(const SimcEvent& evt) const {
    // Port from physics_pion.f
    // Pion electroproduction cross section
    
    // Use MAID model or simple parameterization
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
    // Check pion production kinematics
    
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // W must be above pion production threshold
    double W_threshold = Mp + Mpi;  // ~ 1078 MeV
    if (evt.W < W_threshold) return false;
    
    // W should be reasonable (not too high)
    if (evt.W > 3000.0) return false;  // 3 GeV
    
    return true;
}

double PionCrossSection::HadronicTensor(const SimcEvent& /* evt */) const {
    // Hadronic tensor for pion production
    // Placeholder for full implementation
    
    return 1.0;
}

// ============================================================================
// KaonCrossSection Implementation
// ============================================================================

KaonCrossSection::KaonCrossSection(KaonType type)
    : kaon_type_(type) {
}

double KaonCrossSection::Calculate(const SimcEvent& evt) const {
    // Port from physics_kaon.f
    // Kaon electroproduction cross section
    
    double W_GeV = evt.W / 1.0e3;
    
    // Kaon production threshold
    double W_threshold = (Mp + Mk) / 1.0e3;  // GeV
    
    if (W_GeV < W_threshold) return 0.0;
    
    // Simple Regge model parameterization
    // Full implementation would use Saghai model or MAID
    
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
    // Check kaon production kinematics
    
    if (evt.Ein <= 0.0 || evt.e_E <= 0.0 || evt.p_E <= 0.0) return false;
    if (evt.Q2 <= 0.0) return false;
    
    // W must be above kaon production threshold
    double W_threshold = Mp + Mk;  // ~ 1432 MeV
    if (evt.W < W_threshold) return false;
    
    // W should be reasonable
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
            return std::make_unique<QuasiElasticCrossSection>();
            
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
    // (Could add pion type, kaon type, etc. from config)
    
    return xs;
}

// ============================================================================
// RadiativeCorrections Implementation
// ============================================================================

RadiativeCorrections::RadiativeCorrections(double bt) : bt_(bt) {
}

double RadiativeCorrections::GetCorrectionFactor(double Ein, double Ee, 
                                                 double theta, int Z) const {
    // Port from radc.f
    // Calculate total radiative correction factor
    
    double internal = InternalCorrection(Ein, Ee, theta);
    double external = ExternalCorrection(Ein, Ee, Z, 0.1);  // 0.1 g/cm^2 thickness
    
    // Total correction (multiplicative)
    return (1.0 + internal) * (1.0 + external);
}

double RadiativeCorrections::InternalCorrection(double Ein, double Ee, double theta) const {
    // Port from radc.f - basicrad()
    // Internal bremsstrahlung correction (Schwinger term + vertex correction)
    
    double Q2 = 2.0 * Ein * Ee * (1.0 - std::cos(theta));
    
    // Leading order radiative correction (Schwinger)
    double delta_schwinger = (alpha / pi) * (std::log(Q2 / Me2) - 1.0);
    
    // Vertex correction
    double delta_vertex = -(alpha / pi) * (pi*pi / 6.0);
    
    // Vacuum polarization
    double delta_vacuum = (alpha / (3.0 * pi)) * std::log(Q2 / Me2);
    
    return delta_schwinger + delta_vertex + delta_vacuum;
}

double RadiativeCorrections::ExternalCorrection(double Ein, double Ee, 
                                                int Z, double thickness) const {
    // Port from brem.f
    // External bremsstrahlung in target
    
    // Radiation length for target (approximate)
    double X0 = 36.1 / Z;  // g/cm^2 (rough approximation)
    double t = thickness / X0;  // Thickness in radiation lengths
    
    // Energy loss due to bremsstrahlung
    double E_avg = (Ein + Ee) / 2.0;
    double delta_E = bt_ * E_avg * t;
    
    // Fractional correction
    return delta_E / Ee;
}

double RadiativeCorrections::RadiativeTail(double Ein, double Ee, double theta) const {
    // Radiative tail contribution
    // This is a simplified version - full implementation would use
    // peaking approximation and integrate over radiated photon spectrum
    
    double Q2 = 2.0 * Ein * Ee * (1.0 - std::cos(theta));
    double nu = Ein - Ee;
    
    // Effective radiator thickness
    double b = (alpha / pi) * (std::log(Q2 / Me2) - 1.0);
    
    // Tail integral (simplified)
    double tail = b * nu / Ein;
    
    return tail;
}

} // namespace simc
