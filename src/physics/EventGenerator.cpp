// src/physics/EventGenerator.cpp
// Event generation for SIMC Monte Carlo
// Ported from event.f

#include "simc/EventGenerator.h"
#include "simc/EnergyLoss.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace simc {

// ============================================================================
// Constructor and Initialization
// ============================================================================

EventGenerator::EventGenerator(
    const ConfigManager& config,
    std::shared_ptr<RandomGenerator> random,
    std::shared_ptr<CrossSection> cross_section,
    std::shared_ptr<SpectrometerOptics> optics
)
    : random_(random)
    , cross_section_(cross_section)
    , optics_(optics)
{
    // Load configuration
    auto gen_config = config.GetSection("generation");
    auto target_config = config.GetSection("target");
    auto beam_config = config.GetSection("beam");
    auto spec_e_config = config.GetSection("spectrometer_electron");
    auto spec_p_config = config.GetSection("spectrometer_hadron");
    
    // Beam parameters
    beam_energy_ = beam_config["energy"].get<double>();
    beam_energy_spread_ = beam_config["energy_spread"].get<double>(0.0);
    
    // Generation limits
    gen_limits_.electron.delta_min = gen_config["electron"]["delta_min"].get<double>();
    gen_limits_.electron.delta_max = gen_config["electron"]["delta_max"].get<double>();
    gen_limits_.electron.yptar_min = gen_config["electron"]["yptar_min"].get<double>();
    gen_limits_.electron.yptar_max = gen_config["electron"]["yptar_max"].get<double>();
    gen_limits_.electron.xptar_min = gen_config["electron"]["xptar_min"].get<double>();
    gen_limits_.electron.xptar_max = gen_config["electron"]["xptar_max"].get<double>();
    gen_limits_.electron.E_min = gen_config["electron"]["E_min"].get<double>();
    gen_limits_.electron.E_max = gen_config["electron"]["E_max"].get<double>();
    
    gen_limits_.hadron.delta_min = gen_config["hadron"]["delta_min"].get<double>();
    gen_limits_.hadron.delta_max = gen_config["hadron"]["delta_max"].get<double>();
    gen_limits_.hadron.yptar_min = gen_config["hadron"]["yptar_min"].get<double>();
    gen_limits_.hadron.yptar_max = gen_config["hadron"]["yptar_max"].get<double>();
    gen_limits_.hadron.xptar_min = gen_config["hadron"]["xptar_min"].get<double>();
    gen_limits_.hadron.xptar_max = gen_config["hadron"]["xptar_max"].get<double>();
    gen_limits_.hadron.E_min = gen_config["hadron"]["E_min"].get<double>();
    gen_limits_.hadron.E_max = gen_config["hadron"]["E_max"].get<double>();
    
    // Beam widths
    gen_limits_.xwid = gen_config["beam_xwidth"].get<double>(0.01); // cm, 1 sigma
    gen_limits_.ywid = gen_config["beam_ywidth"].get<double>(0.01); // cm, 1 sigma
    
    // Target properties
    target_props_.mass = target_config["mass"].get<double>();
    target_props_.Z = target_config["Z"].get<int>();
    target_props_.A = target_config["A"].get<int>();
    target_props_.length = target_config["length"].get<double>();
    target_props_.x_offset = target_config["x_offset"].get<double>(0.0);
    target_props_.y_offset = target_config["y_offset"].get<double>(0.0);
    target_props_.z_offset = target_config["z_offset"].get<double>(0.0);
    target_props_.angle = target_config["angle"].get<double>(0.0);
    target_props_.raster_pattern = target_config["raster_pattern"].get<int>(0);
    target_props_.raster_x = target_config["raster_x"].get<double>(0.0);
    target_props_.raster_y = target_config["raster_y"].get<double>(0.0);
    
    // Spectrometer settings
    spec_electron_.P = spec_e_config["momentum"].get<double>();
    spec_electron_.theta = spec_e_config["angle"].get<double>() * kDegToRad;
    spec_electron_.phi = spec_e_config["phi"].get<double>(0.0) * kDegToRad;
    
    spec_hadron_.P = spec_p_config["momentum"].get<double>();
    spec_hadron_.theta = spec_p_config["angle"].get<double>() * kDegToRad;
    spec_hadron_.phi = spec_p_config["phi"].get<double>(0.0) * kDegToRad;
    
    // Flags
    use_energy_loss_ = gen_config["use_energy_loss"].get<bool>(true);
    use_coulomb_ = gen_config["use_coulomb"].get<bool>(true);
    use_radiative_ = gen_config["use_radiative"].get<bool>(false);
    
    // Reaction type
    std::string reaction = gen_config["reaction_type"].get<std::string>("H(e,e'p)");
    if (reaction == "H(e,e'p)") {
        reaction_type_ = ReactionType::kHydElastic;
    } else if (reaction == "D(e,e'p)") {
        reaction_type_ = ReactionType::kDeuterium;
    } else if (reaction == "A(e,e'p)") {
        reaction_type_ = ReactionType::kHeavy;
    } else if (reaction == "pion") {
        reaction_type_ = ReactionType::kPion;
    } else if (reaction == "kaon") {
        reaction_type_ = ReactionType::kKaon;
    } else {
        throw std::runtime_error("Unknown reaction type: " + reaction);
    }
}

bool EventGenerator::Initialize() {
    if (initialized_) {
        return true;
    }
    
    // Calculate average beam energy at vertex including energy loss
    // From event.f init.f calculation of Ebeam_vertex_ave
    double eloss_ave = 0.0;
    if (use_energy_loss_) {
        // Use middle of target
        double teff;
        TripThruTarget(1, 0.0, beam_energy_, 0.0, eloss_ave, teff);
    }
    
    double coulomb_ave = 0.0;
    if (use_coulomb_) {
        // Average Coulomb correction (from target.f)
        coulomb_ave = target_props_.coulomb_constant * 
                      (3.0 - std::pow(0.5, 2.0/3.0));
    }
    
    beam_energy_vertex_ = beam_energy_ + coulomb_ave - eloss_ave;
    
    std::cout << "EventGenerator initialized:" << std::endl;
    std::cout << "  Beam energy: " << beam_energy_ << " MeV" << std::endl;
    std::cout << "  Beam energy at vertex: " << beam_energy_vertex_ << " MeV" << std::endl;
    std::cout << "  Reaction type: " << static_cast<int>(reaction_type_) << std::endl;
    std::cout << "  Target: Z=" << target_props_.Z 
              << " A=" << target_props_.A 
              << " M=" << target_props_.mass << " MeV" << std::endl;
    
    initialized_ = true;
    return true;
}

// ============================================================================
// Main Event Generation
// ============================================================================

bool EventGenerator::GenerateEvent(SimcEvent& vertex, MainEvent& main) {
    // From event.f generate() subroutine
    
    if (!initialized_) {
        throw std::runtime_error("EventGenerator not initialized");
    }
    
    n_generated_++;
    
    // Initialize
    vertex.Clear();
    main.Clear();
    bool success = false;
    main.gen_weight = 1.0;
    
    // Generate vertex position
    GenerateVertex(main.target);
    
    // Calculate energy loss to vertex and Coulomb correction
    // From event.f lines 230-250
    double eloss_to_vertex, teff;
    TripThruTarget(1, main.target.z - target_props_.z_offset, 
                   beam_energy_, 0.0, eloss_to_vertex, teff);
    
    if (!use_energy_loss_) {
        eloss_to_vertex = 0.0;
    }
    
    main.target.Eloss[0] = eloss_to_vertex;
    main.target.teff[0] = teff;
    
    // Coulomb correction
    if (use_coulomb_) {
        main.target.Coulomb = target_props_.coulomb_constant;
    } else {
        main.target.Coulomb = 0.0;
    }
    
    // Beam energy at vertex with fluctuations
    // From event.f lines 255-260
    vertex.Ein = beam_energy_ + 
                 (random_->Uniform() - 0.5) * beam_energy_spread_ +
                 main.target.Coulomb - main.target.Eloss[0];
    
    // Energy shifts for comparison with limits
    main.Ein_shift = vertex.Ein - beam_energy_vertex_;
    main.Ee_shift = main.target.Coulomb - target_props_.coulomb_ave;
    
    // Generate phase space
    if (!GeneratePhaseSpace(vertex, main)) {
        return false;
    }
    
    // Complete event kinematics
    if (!CompleteEvent(vertex, main)) {
        return false;
    }
    
    // Calculate energy loss for outgoing particles
    // From event.f lines 890-905
    TripThruTarget(2, main.target.z - target_props_.z_offset,
                   vertex.e_E, vertex.e_theta, 
                   main.target.Eloss[1], main.target.teff[1]);
    
    TripThruTarget(3, main.target.z - target_props_.z_offset,
                   vertex.p_E, vertex.p_theta,
                   main.target.Eloss[2], main.target.teff[2]);
    
    if (!use_energy_loss_) {
        main.target.Eloss[1] = 0.0;
        main.target.Eloss[2] = 0.0;
    }
    
    // Check if event passes cuts
    success = PassesGenerationCuts(vertex);
    
    if (success) {
        n_accepted_++;
    }
    
    main.success = success;
    return success;
}

// ============================================================================
// Vertex Generation
// ============================================================================

void EventGenerator::GenerateVertex(TargetInfo& target) {
    // From event.f generate() lines 150-230
    
    // Generate beam position with Gaussian profile (±3σ)
    // From event.f: main.target.x = gauss1(nsig_max)*gen.xwid+targ.xoffset
    target.x = random_->Gauss(0.0, gen_limits_.xwid, kNSigmaMax) + 
               target_props_.x_offset;
    target.y = random_->Gauss(0.0, gen_limits_.ywid, kNSigmaMax) + 
               target_props_.y_offset;
    
    // Add raster pattern
    // From event.f lines 180-210
    double raster_x = 0.0;
    double raster_y = 0.0;
    
    if (target_props_.raster_pattern == 1) {
        // Old bedpost raster - square with cosine distribution
        double t3 = random_->Uniform() * kPi;
        double t4 = random_->Uniform() * kPi;
        raster_x = std::cos(t3) * target_props_.raster_x;
        raster_y = std::cos(t4) * target_props_.raster_y;
        
    } else if (target_props_.raster_pattern == 2) {
        // Circular raster - uniform in annulus
        double t3 = random_->Uniform() * 2.0 * kPi;
        double t4 = std::sqrt(random_->Uniform()) * 
                    (target_props_.raster_y - target_props_.raster_x) + 
                    target_props_.raster_x;
        raster_x = std::cos(t3) * t4;
        raster_y = std::sin(t3) * t4;
        
    } else if (target_props_.raster_pattern == 3) {
        // New flat square raster - uniform
        double t3 = 2.0 * random_->Uniform() - 1.0;
        double t4 = 2.0 * random_->Uniform() - 1.0;
        raster_x = target_props_.raster_x * t3;
        raster_y = target_props_.raster_y * t4;
    }
    
    target.x += raster_x;
    target.y += raster_y;
    target.raster_x = raster_x;
    target.raster_y = raster_y;
    
    // Position along target (uniform)
    // From event.f: main.target.z = (0.5-grnd())*targ.length+targ.zoffset
    target.z = (0.5 - random_->Uniform()) * target_props_.length + 
               target_props_.z_offset;
}

// ============================================================================
// Phase Space Generation
// ============================================================================

bool EventGenerator::GeneratePhaseSpace(SimcEvent& vertex, MainEvent& main) {
    // From event.f generate() lines 260-400
    
    // Always generate electron angles
    // From event.f lines 260-265
    vertex.e_yptar = gen_limits_.electron.yptar_min + 
                     random_->Uniform() * (gen_limits_.electron.yptar_max - 
                                          gen_limits_.electron.yptar_min);
    vertex.e_xptar = gen_limits_.electron.xptar_min + 
                     random_->Uniform() * (gen_limits_.electron.xptar_max - 
                                          gen_limits_.electron.xptar_min);
    
    // Generate hadron angles for all except H(e,e'p)
    // From event.f lines 267-275
    if (reaction_type_ != ReactionType::kHydElastic) {
        vertex.p_yptar = gen_limits_.hadron.yptar_min + 
                         random_->Uniform() * (gen_limits_.hadron.yptar_max - 
                                              gen_limits_.hadron.yptar_min);
        vertex.p_xptar = gen_limits_.hadron.xptar_min + 
                         random_->Uniform() * (gen_limits_.hadron.xptar_max - 
                                              gen_limits_.hadron.xptar_min);
    }
    
    // Generate energies based on reaction type
    if (reaction_type_ == ReactionType::kDeuterium ||
        reaction_type_ == ReactionType::kHeavy) {
        
        // Generate electron energy
        // From event.f lines 290-310
        double Emin = gen_limits_.electron.E_min;
        double Emax = gen_limits_.electron.E_max;
        
        // Apply sumEgen constraint if needed
        if (reaction_type_ == ReactionType::kHeavy) {
            Emin = std::max(Emin, gen_limits_.sumEgen_min - gen_limits_.hadron.E_max);
            Emax = std::min(Emax, gen_limits_.sumEgen_max - gen_limits_.hadron.E_min);
        }
        
        if (Emin > Emax) {
            return false;
        }
        
        // Calculate generation weight
        main.gen_weight *= (Emax - Emin) / 
                          (gen_limits_.electron.E_max - gen_limits_.electron.E_min);
        
        vertex.e_E = Emin + random_->Uniform() * (Emax - Emin);
        vertex.e_P = vertex.e_E; // Ultra-relativistic electron
        vertex.e_delta = 100.0 * (vertex.e_P - spec_electron_.P) / spec_electron_.P;
    }
    
    // For H(e,e'p), electron energy will be calculated in SolveHydrogenElastic()
    
    // Convert spectrometer angles to physics angles
    // From event.f lines 320-330
    PhysicsAngles(spec_electron_.theta, spec_electron_.phi,
                  vertex.e_xptar, vertex.e_yptar,
                  vertex.e_theta, vertex.e_phi);
    
    if (reaction_type_ != ReactionType::kHydElastic) {
        PhysicsAngles(spec_hadron_.theta, spec_hadron_.phi,
                      vertex.p_xptar, vertex.p_yptar,
                      vertex.p_theta, vertex.p_phi);
    }
    
    return true;
}

// ============================================================================
// Complete Event Kinematics
// ============================================================================

bool EventGenerator::CompleteEvent(SimcEvent& vertex, MainEvent& main) {
    // From event.f complete_ev() subroutine
    
    main.jacobian = 1.0;
    
    // Calculate unit vectors
    // From event.f complete_ev() lines 495-510
    vertex.ue_x = std::sin(vertex.e_theta) * std::cos(vertex.e_phi);
    vertex.ue_y = std::sin(vertex.e_theta) * std::sin(vertex.e_phi);
    vertex.ue_z = std::cos(vertex.e_theta);
    
    if (reaction_type_ != ReactionType::kHydElastic) {
        vertex.up_x = std::sin(vertex.p_theta) * std::cos(vertex.p_phi);
        vertex.up_y = std::sin(vertex.p_theta) * std::sin(vertex.p_phi);
        vertex.up_z = std::cos(vertex.p_theta);
    }
    
    // Reaction-specific kinematics
    bool success = false;
    
    switch (reaction_type_) {
        case ReactionType::kHydElastic:
            success = SolveHydrogenElastic(vertex);
            break;
            
        case ReactionType::kDeuterium:
            success = SolveDeuterium(vertex, main);
            break;
            
        // Other reactions will be added in later iterations
        default:
            throw std::runtime_error("Reaction type not yet implemented");
    }
    
    if (!success) {
        return false;
    }
    
    // Calculate basic kinematics
    CalculateBasicKinematics(vertex);
    
    // Calculate missing momentum
    CalculateMissingMomentum(vertex);
    
    // Calculate jacobian
    main.jacobian = CalculateJacobian(vertex, main);
    
    // Validate physics
    return ValidatePhysics(vertex);
}

// ============================================================================
// Basic Kinematics Calculations
// ============================================================================

void EventGenerator::CalculateBasicKinematics(SimcEvent& event) const {
    // From event.f complete_ev() lines 500-520
    
    // Energy transfer
    event.nu = event.Ein - event.e_E;
    
    // Four-momentum transfer squared
    event.Q2 = 2.0 * event.Ein * event.e_E * (1.0 - event.ue_z);
    
    // Three-momentum transfer
    event.q = std::sqrt(event.Q2 + event.nu * event.nu);
    
    // Bjorken x
    event.xbj = event.Q2 / (2.0 * kProtonMass * event.nu);
    
    // q unit vector
    event.uq_x = -event.e_P * event.ue_x / event.q;
    event.uq_y = -event.e_P * event.ue_y / event.q;
    event.uq_z = (event.Ein - event.e_P * event.ue_z) / event.q;
    
    // Verify normalization
    double uq_norm = event.uq_x * event.uq_x + 
                     event.uq_y * event.uq_y + 
                     event.uq_z * event.uq_z;
    if (std::abs(uq_norm - 1.0) > 0.01) {
        std::cerr << "Warning: q vector normalization error: " << uq_norm << std::endl;
    }
}

void EventGenerator::CalculateMissingMomentum(SimcEvent& event) const {
    // From event.f complete_ev() lines 750-800
    
    // Missing momentum vector: P_miss = P_hadron - q
    event.Pmx = event.p_P * event.up_x - event.q * event.uq_x;
    event.Pmy = event.p_P * event.up_y - event.q * event.uq_y;
    event.Pmz = event.p_P * event.up_z - event.q * event.uq_z;
    event.Pmiss = std::sqrt(event.Pmx * event.Pmx + 
                           event.Pmy * event.Pmy + 
                           event.Pmz * event.Pmz);
    
    // Missing energy
    event.Emiss = event.nu + target_props_.mass - event.p_E;
    
    // Out-of-plane unit vector: oop = q_hat × z_hat (but normalized)
    double oop_x = -event.uq_y;  // = q_y * (z × y) = -q_y * x
    double oop_y =  event.uq_x;  // = q_x * (z × x) = q_x * y
    double oop_norm = std::sqrt(oop_x * oop_x + oop_y * oop_y);
    
    // Parallel component (along q)
    event.PmPar = event.Pmx * event.uq_x + 
                  event.Pmy * event.uq_y + 
                  event.Pmz * event.uq_z;
    
    // Out-of-plane component
    if (oop_norm > 0.0) {
        event.PmOop = (event.Pmx * oop_x + event.Pmy * oop_y) / oop_norm;
    } else {
        event.PmOop = 0.0;
    }
    
    // Perpendicular component (what's left)
    double pm_sq = event.Pmiss * event.Pmiss - 
                   event.PmPar * event.PmPar - 
                   event.PmOop * event.PmOop;
    event.PmPer = std::sqrt(std::max(0.0, pm_sq));
}


// ============================================================================
// Hydrogen Elastic Kinematics
// ============================================================================

bool EventGenerator::SolveHydrogenElastic(SimcEvent& event) const {
    // From event.f complete_ev() lines 525-540
    
    // For elastic scattering: E_e' = E_in * M / (M + E_in * (1 - cos(theta)))
    double denominator = kProtonMass + event.Ein * (1.0 - event.ue_z);
    event.e_E = event.Ein * kProtonMass / denominator;
    
    // Check for unphysical solution
    if (event.e_E > event.Ein) {
        return false;
    }
    
    event.e_P = event.e_E; // Ultra-relativistic
    event.e_delta = 100.0 * (event.e_P - spec_electron_.P) / spec_electron_.P;
    
    // For elastic scattering, proton momentum equals q
    event.Em = 0.0;
    event.Pm = 0.0;
    event.Mrec = 0.0;
    event.Trec = 0.0;
    
    // Proton direction = q direction
    event.up_x = event.uq_x;
    event.up_y = event.uq_y;
    event.up_z = event.uq_z;
    
    // Proton momentum
    event.p_P = event.q;
    event.p_E = std::sqrt(event.p_P * event.p_P + kProtonMass * kProtonMass);
    event.p_delta = 100.0 * (event.p_P - spec_hadron_.P) / spec_hadron_.P;
    
    // Calculate physics angles for proton
    event.p_theta = std::acos(event.up_z);
    
    double sin_theta = std::sin(event.p_theta);
    if (sin_theta > 0.0) {
        double cos_phi = event.up_x / sin_theta;
        if (std::abs(cos_phi) > 1.0) {
            // Numerical issue
            cos_phi = std::copysign(1.0, cos_phi);
        }
        event.p_phi = std::atan2(event.up_y, event.up_x);
        if (event.p_phi < 0.0) {
            event.p_phi += 2.0 * kPi;
        }
    } else {
        event.p_phi = 0.0;
    }
    
    // Convert back to spectrometer angles
    SpectrometerAngles(spec_hadron_.theta, spec_hadron_.phi,
                       event.p_theta, event.p_phi,
                       event.p_xptar, event.p_yptar);
    
    return true;
}

// ============================================================================
// Deuterium Kinematics with Jacobian
// ============================================================================

bool EventGenerator::SolveDeuterium(SimcEvent& event, MainEvent& main) const {
    // From event.f complete_ev() lines 545-590
    
    // Fixed missing energy (binding energy of deuteron)
    event.Em = kProtonMass + kNeutronMass - target_props_.mass; // ~2.2249 MeV
    event.Mrec = target_props_.mass - kProtonMass + event.Em;   // = neutron mass
    
    // Solve quadratic equation for proton energy
    // From momentum and energy conservation:
    // (E_d + nu - E_p)² = M_n² + (q - P_p)²
    
    // Calculate coefficients
    double a = -event.q * (event.uq_x * event.up_x + 
                          event.uq_y * event.up_y + 
                          event.uq_z * event.up_z);
    double b = event.q * event.q;
    double c = event.nu + target_props_.mass;
    double t = c * c - b + kProtonMass * kProtonMass - event.Mrec * event.Mrec;
    
    double QA = 4.0 * (a * a - c * c);
    double QB = 4.0 * c * t;
    double QC = -4.0 * a * a * kProtonMass * kProtonMass - t * t;
    
    double radical = QB * QB - 4.0 * QA * QC;
    
    if (radical < 0.0) {
        return false; // No physical solution
    }
    
    // Take the higher solution (forward-going proton)
    event.p_E = (-QB - std::sqrt(radical)) / (2.0 * QA);
    
    // Check for unphysical solution
    if (event.p_E <= kProtonMass) {
        return false;
    }
    
    // Check second solution (should be outside acceptance)
    double E_proton_2 = (-QB + std::sqrt(radical)) / (2.0 * QA);
    if (E_proton_2 > gen_limits_.hadron.E_min && 
        E_proton_2 < gen_limits_.hadron.E_max) {
        // This is problematic - both solutions in acceptance
        std::cerr << "Warning: Both quadratic solutions in acceptance: " 
                  << event.p_E << ", " << E_proton_2 << std::endl;
    }
    
    // Calculate proton momentum and delta
    event.p_P = std::sqrt(event.p_E * event.p_E - kProtonMass * kProtonMass);
    event.p_delta = 100.0 * (event.p_P - spec_hadron_.P) / spec_hadron_.P;
    
    // Calculate Jacobian |dE_p/dE_m|
    // From event.f lines 580-585
    double jac_num = t * (c - event.p_E) + 2.0 * c * event.p_E * (event.p_E - c);
    double jac_den = 2.0 * (a * a - c * c) * event.p_E + c * t;
    main.jacobian *= std::abs(jac_num / jac_den);
    
    // Missing momentum
    event.Pm = event.Pmiss; // Will be calculated in CalculateMissingMomentum
    
    // Recoil kinetic energy
    event.Trec = std::sqrt(event.Mrec * event.Mrec + event.Pm * event.Pm) - 
                 event.Mrec;
    
    return true;
}

// ============================================================================
// Angle Conversions
// ============================================================================

void EventGenerator::PhysicsAngles(
    Double_t theta0, Double_t phi0,
    Double_t xptar, Double_t yptar,
    Double_t& theta, Double_t& phi
) const {
    // From event.f physics_angles() subroutine
    
    double cos_theta0 = std::cos(theta0);
    double sin_theta0 = std::sin(theta0);
    double sin_phi0 = std::sin(phi0);
    double cos_phi0 = std::cos(phi0);
    
    double r = std::sqrt(1.0 + xptar * xptar + yptar * yptar);
    
    // Warn if phi0 is not ±π/2 (formulas assume this)
    if (std::abs(cos_phi0) > 0.0001) {
        std::cerr << "Warning: phi0 != ±π/2 may give incorrect results. phi0 = " 
                  << phi0 * kRadToDeg << " degrees" << std::endl;
    }
    
    // Calculate theta
    double cos_theta_arg = (cos_theta0 - yptar * sin_theta0 * sin_phi0) / r;
    
    // Clamp to valid range for acos
    if (cos_theta_arg > 1.0) cos_theta_arg = 1.0;
    if (cos_theta_arg < -1.0) cos_theta_arg = -1.0;
    
    theta = std::acos(cos_theta_arg);
    
    // Calculate phi
    if (xptar != 0.0) {
        double phi_num = yptar * cos_theta0 + sin_theta0 * sin_phi0;
        phi = std::atan(phi_num / xptar);
        
        // Adjust to 0-180 degrees
        if (phi <= 0.0) {
            phi += kPi;
        }
        
        // Add π for HMS (sin_phi0 < 0)
        if (sin_phi0 < 0.0) {
            phi += kPi;
        }
    } else {
        phi = phi0;
    }
}

void EventGenerator::SpectrometerAngles(
    Double_t theta0, Double_t phi0,
    Double_t theta, Double_t phi,
    Double_t& xptar, Double_t& yptar
) const {
    // From event.f spectrometer_angles() subroutine
    
    // Unit vector in physics frame
    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);
    
    // Central spectrometer direction
    double x0 = std::sin(theta0) * std::cos(phi0);
    double y0 = std::sin(theta0) * std::sin(phi0);
    double z0 = std::cos(theta0);
    
    // Cos of angle between particle and spectrometer central ray
    double cos_dtheta = x * x0 + y * y0 + z * z0;
    
    // Project to spectrometer plane
    xptar = x / cos_dtheta;
    
    double yptar_mag = std::sqrt(1.0 / (cos_dtheta * cos_dtheta) - 1.0 - xptar * xptar);
    
    // Determine sign of yptar
    double y_event = y / cos_dtheta;
    if (y_event < y0) {
        yptar = -yptar_mag;
    } else {
        yptar = yptar_mag;
    }
}

// ============================================================================
// Jacobian Calculations
// ============================================================================

Double_t EventGenerator::CalculateJacobian(
    const SimcEvent& vertex, 
    const MainEvent& main
) const {
    // From event.f complete_ev() lines 850-870
    
    double jacobian = main.jacobian; // May already include deuterium jacobian
    
    // Angle jacobian for electron: 1/cos³(θ-θ₀)
    double r_e = std::sqrt(1.0 + vertex.e_yptar * vertex.e_yptar + 
                                  vertex.e_xptar * vertex.e_xptar);
    jacobian /= (r_e * r_e * r_e);
    
    // Angle jacobian for hadron (for all except phase space)
    if (reaction_type_ != ReactionType::kPhaseSpace) {
        if (reaction_type_ != ReactionType::kHydElastic || 
            reaction_type_ == ReactionType::kDeuterium ||
            reaction_type_ == ReactionType::kHeavy) {
            
            double r_p = std::sqrt(1.0 + vertex.p_yptar * vertex.p_yptar + 
                                        vertex.p_xptar * vertex.p_xptar);
            jacobian /= (r_p * r_p * r_p);
        }
    }
    
    return jacobian;
}

// ============================================================================
// Helper Methods
// ============================================================================

void EventGenerator::TripThruTarget(
    Int_t particle_id,
    Double_t z_position,
    Double_t energy,
    Double_t theta,
    Double_t& eloss,
    Double_t& teff
) const {
    // Calculate energy loss through target
    // From event.f trip_thru_target() (in target.f)
    
    // Path length through target accounting for angle
    double path_length = std::abs(z_position) / std::cos(theta);
    
    // Effective thickness (g/cm²)
    teff = target_props_.density * path_length;
    
    // Use EnergyLoss class for calculation
    double particle_mass;
    if (particle_id == 1 || particle_id == 2) {
        particle_mass = kElectronMass;
    } else {
        particle_mass = kProtonMass; // Or use Mh for hadron
    }
    
    // Simple energy loss calculation (will use EnergyLoss class in full version)
    // For now, use dE/dx * thickness
    eloss = EnergyLoss::BetheBloch(energy, particle_mass, 
                                   target_props_.Z, target_props_.A) * teff;
}

Vector3D EventGenerator::UnitVector(Double_t theta, Double_t phi) const {
    return {
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
        std::cos(theta)
    };
}

bool EventGenerator::PassesGenerationCuts(const SimcEvent& event) const {
    // Check if event is within generation limits
    
    // Electron
    if (event.e_delta < gen_limits_.electron.delta_min || 
        event.e_delta > gen_limits_.electron.delta_max) {
        return false;
    }
    
    if (event.e_xptar < gen_limits_.electron.xptar_min || 
        event.e_xptar > gen_limits_.electron.xptar_max) {
        return false;
    }
    
    if (event.e_yptar < gen_limits_.electron.yptar_min || 
        event.e_yptar > gen_limits_.electron.yptar_max) {
        return false;
    }
    
    // Hadron (for non-elastic)
    if (reaction_type_ != ReactionType::kHydElastic) {
        if (event.p_delta < gen_limits_.hadron.delta_min || 
            event.p_delta > gen_limits_.hadron.delta_max) {
            return false;
        }
        
        if (event.p_xptar < gen_limits_.hadron.xptar_min || 
            event.p_xptar > gen_limits_.hadron.xptar_max) {
            return false;
        }
        
        if (event.p_yptar < gen_limits_.hadron.yptar_min || 
            event.p_yptar > gen_limits_.hadron.yptar_max) {
            return false;
        }
    }
    
    return true;
}

bool EventGenerator::ValidatePhysics(const SimcEvent& event) const {
    // Check for unphysical values
    
    // Energy conservation
    if (event.e_E > event.Ein) {
        return false;
    }
    
    // Q² should be positive
    if (event.Q2 < 0.0) {
        return false;
    }
    
    // Particle energies must be above mass
    if (event.e_E < kElectronMass) {
        return false;
    }
    
    if (event.p_E < kProtonMass) {
        return false;
    }
    
    return true;
}

void EventGenerator::ResetStatistics() {
    n_generated_ = 0;
    n_accepted_ = 0;
}

} // namespace simc
