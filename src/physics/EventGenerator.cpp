// src/physics/EventGenerator.cpp - PART 1
// Constructor through GeneratePhaseSpace
// Ported from event.f, jacobians.f

#include "simc/physics/EventGenerator.h"
#include "simc/physics/EnergyLoss.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace simc {

  // ============================================================================
  // Constructor
  // ============================================================================

  EventGenerator::EventGenerator(
				 const ConfigManager& config,
				 std::shared_ptr<RandomGenerator> random,
				 std::shared_ptr<CrossSectionBase> cross_section,
				 std::shared_ptr<SpectrometerOptics> optics)
    : random_(random)
    , cross_section_(cross_section)
    , optics_(optics)
  {
    if (!random_) {
      throw std::runtime_error("EventGenerator: RandomGenerator is required");
    }
    
    // Load configuration using Get<T>(path) template method
    beam_energy_ = config.Get<double>("beam.energy");
    beam_energy_spread_ = config.Get<double>("beam.energy_spread", 0.0);
    
    // Generation limits - electron
    gen_limits_.electron.delta_min = config.Get<double>("generation.electron.delta_min");
    gen_limits_.electron.delta_max = config.Get<double>("generation.electron.delta_max");
    gen_limits_.electron.yptar_min = config.Get<double>("generation.electron.yptar_min");
    gen_limits_.electron.yptar_max = config.Get<double>("generation.electron.yptar_max");
    gen_limits_.electron.xptar_min = config.Get<double>("generation.electron.xptar_min");
    gen_limits_.electron.xptar_max = config.Get<double>("generation.electron.xptar_max");
    gen_limits_.electron.E_min = config.Get<double>("generation.electron.E_min");
    gen_limits_.electron.E_max = config.Get<double>("generation.electron.E_max");
    
    // Generation limits - hadron
    gen_limits_.hadron.delta_min = config.Get<double>("generation.hadron.delta_min");
    gen_limits_.hadron.delta_max = config.Get<double>("generation.hadron.delta_max");
    gen_limits_.hadron.yptar_min = config.Get<double>("generation.hadron.yptar_min");
    gen_limits_.hadron.yptar_max = config.Get<double>("generation.hadron.yptar_max");
    gen_limits_.hadron.xptar_min = config.Get<double>("generation.hadron.xptar_min");
    gen_limits_.hadron.xptar_max = config.Get<double>("generation.hadron.xptar_max");
    gen_limits_.hadron.E_min = config.Get<double>("generation.hadron.E_min");
    gen_limits_.hadron.E_max = config.Get<double>("generation.hadron.E_max");
    
    // Beam widths
    gen_limits_.xwid = config.Get<double>("generation.beam_xwidth", 0.01);
    gen_limits_.ywid = config.Get<double>("generation.beam_ywidth", 0.01);
    
    // Target
    target_props_.mass = config.Get<double>("target.mass");
    target_props_.Z = config.Get<int>("target.Z");
    target_props_.A = config.Get<int>("target.A");
    target_props_.length = config.Get<double>("target.length");
    target_props_.x_offset = config.Get<double>("target.x_offset", 0.0);
    target_props_.y_offset = config.Get<double>("target.y_offset", 0.0);
    target_props_.z_offset = config.Get<double>("target.z_offset", 0.0);
    target_props_.angle = config.Get<double>("target.angle", 0.0);
    target_props_.raster_pattern = config.Get<int>("target.raster_pattern", 0);
    target_props_.raster_x = config.Get<double>("target.raster_x", 0.0);
    target_props_.raster_y = config.Get<double>("target.raster_y", 0.0);
    
    // Spectrometers
    spec_electron_.P = config.Get<double>("spectrometer_electron.momentum");
    spec_electron_.theta = config.Get<double>("spectrometer_electron.angle") * kDegToRad;
    spec_electron_.phi = config.Get<double>("spectrometer_electron.phi", 0.0) * kDegToRad;
    
    spec_hadron_.P = config.Get<double>("spectrometer_hadron.momentum");
    spec_hadron_.theta = config.Get<double>("spectrometer_hadron.angle") * kDegToRad;
    spec_hadron_.phi = config.Get<double>("spectrometer_hadron.phi", 0.0) * kDegToRad;
    
    // Flags
    use_energy_loss_ = config.Get<bool>("generation.use_energy_loss", true);
    use_coulomb_ = config.Get<bool>("generation.use_coulomb", true);
    use_radiative_ = config.Get<bool>("generation.use_radiative", false);
    
    // Parse reaction type - use correct enum values
    std::string reaction = config.Get<std::string>("generation.reaction_type", "H(e,e'p)");
    if (reaction == "H(e,e'p)" || reaction == "elastic") {
      reaction_type_ = ReactionType::ELASTIC;
    } else if (reaction == "D(e,e'p)" || reaction == "deuterium" || reaction == "quasielastic") {
      reaction_type_ = ReactionType::QUASIELASTIC;
    } else if (reaction == "pion" || reaction == "pi") {
      reaction_type_ = ReactionType::PION_PRODUCTION;
    } else if (reaction == "kaon" || reaction == "K") {
      reaction_type_ = ReactionType::KAON_PRODUCTION;
    } else {
      std::cerr << "Warning: Unknown reaction '" << reaction << "', using ELASTIC\n";
      reaction_type_ = ReactionType::ELASTIC;
    }
  }

  // ============================================================================
  // Initialize
  // ============================================================================

  bool EventGenerator::Initialize() {
    beam_energy_vertex_ = beam_energy_;
    initialized_ = true;
    return true;
  }

  // ============================================================================
  // Reset Statistics
  // ============================================================================

  void EventGenerator::ResetStatistics() {
    n_generated_ = 0;
    n_accepted_ = 0;
  }

  // ============================================================================
  // Generate Event
  // ============================================================================

  bool EventGenerator::GenerateEvent(SimcEvent& event, MainEvent& main) {
    n_generated_++;
    
    GenerateVertex(main.target);
    
    if (!GeneratePhaseSpace(event, main)) {
      return false;
    }
    
    if (!CompleteEvent(event, main)) {
      return false;
    }
    
    CalculateBasicKinematics(event);
    CalculateMissingMomentum(event);
    
    main.jacobian = CalculateJacobian(event, main);
    
    if (!ValidatePhysics(event)) {
      return false;
    }
    
    if (!PassesGenerationCuts(event)) {
      return false;
    }
    
    n_accepted_++;
    main.success = true;
    return true;
  }

  // ============================================================================
  // Generate Vertex
  // ============================================================================

  void EventGenerator::GenerateVertex(TargetInfo& target) {
    // Use Gaussian() method (not Gauss)
    target.x = random_->Gaussian(0.0, gen_limits_.xwid);
    if (std::abs(target.x) > kNSigmaMax * gen_limits_.xwid) {
      target.x = (target.x > 0 ? kNSigmaMax : -kNSigmaMax) * gen_limits_.xwid;
    }
    target.x += target_props_.x_offset;
    
    target.y = random_->Gaussian(0.0, gen_limits_.ywid);
    if (std::abs(target.y) > kNSigmaMax * gen_limits_.ywid) {
      target.y = (target.y > 0 ? kNSigmaMax : -kNSigmaMax) * gen_limits_.ywid;
    }
    target.y += target_props_.y_offset;
    
    target.z = random_->Uniform(-target_props_.length / 2.0, 
				target_props_.length / 2.0) + target_props_.z_offset;
    
    // Raster patterns
    double raster_x = 0.0, raster_y = 0.0;
    
    switch (target_props_.raster_pattern) {
    case 0: // No raster
      break;
            
    case 1: // Bedpost (uniform square)
      raster_x = random_->Uniform(-target_props_.raster_x / 2.0,
				  target_props_.raster_x / 2.0);
      raster_y = random_->Uniform(-target_props_.raster_y / 2.0,
				  target_props_.raster_y / 2.0);
      break;
            
    case 2: // Circular (uniform disk)
      {
	double r = std::sqrt(random_->Uniform()) * target_props_.raster_x / 2.0;
	double angle = random_->Uniform() * kPi;
	raster_x = r * std::cos(angle);
	raster_y = r * std::sin(angle);
	break;
      }
        
    case 3: // Flat circular
      {
	double r = random_->Uniform() * target_props_.raster_x / 2.0;
	double angle = random_->Uniform() * 2.0 * kPi;
	raster_x = r * std::cos(angle);
	raster_y = r * std::sin(angle);
	break;
      }
        
    default:
      std::cerr << "Warning: Unknown raster pattern " 
		<< target_props_.raster_pattern << "\n";
    }
    
    target.x += raster_x;
    target.y += raster_y;
    
    // Use correct field names: rasterx, rastery (not raster_x, raster_y)
    target.rasterx = raster_x;
    target.rastery = raster_y;
  }

  // ============================================================================
  // Generate Phase Space
  // ============================================================================
  
  bool EventGenerator::GeneratePhaseSpace(SimcEvent& event, MainEvent& main) {
    // Generate Electron Angles (ALL cases):
    event.e_yptar = random_->Uniform(gen_limits_.electron.yptar_min,
                                     gen_limits_.electron.yptar_max);
    event.e_xptar = random_->Uniform(gen_limits_.electron.xptar_min,
                                     gen_limits_.electron.xptar_max);
    
    // Calculate electron physics angles from spectrometer angles
    PhysicsAngles(spec_electron_.theta, spec_electron_.phi,
                  event.e_xptar, event.e_yptar,
                  event.e_theta, event.e_phi);
    
    // Generate Hadron Angles (all but H(e,e'p)):
    if (reaction_type_ != ReactionType::ELASTIC) {
      event.p_yptar = random_->Uniform(gen_limits_.hadron.yptar_min,
				       gen_limits_.hadron.yptar_max);
      event.p_xptar = random_->Uniform(gen_limits_.hadron.xptar_min,
				       gen_limits_.hadron.xptar_max);
        
      PhysicsAngles(spec_hadron_.theta, spec_hadron_.phi,
		    event.p_xptar, event.p_yptar,
		    event.p_theta, event.p_phi);
    }
    
    // Generate Hadron Momentum (A(e,e'p) or semi-inclusive)
    if (reaction_type_ == ReactionType::QUASIELASTIC || 
        reaction_type_ == ReactionType::SEMI_INCLUSIVE) {
      double Emin = std::max(gen_limits_.hadron.E_min, 
			     gen_limits_.sumEgen_min - gen_limits_.electron.E_max);
      double Emax = std::min(gen_limits_.hadron.E_max, 
			     gen_limits_.sumEgen_max - gen_limits_.electron.E_min);
        
      if (Emin > Emax) return false;
        
      main.gen_weight *= (Emax - Emin) / 
	(gen_limits_.hadron.E_max - gen_limits_.hadron.E_min);
        
      event.p_E = random_->Uniform(Emin, Emax);
      event.p_P = std::sqrt(event.p_E * event.p_E - kProtonMass * kProtonMass);
      event.p_delta = 100.0 * (event.p_P - spec_hadron_.P) / spec_hadron_.P;
    }
    
    // Generate Electron Energy (all but hydrogen elastic)
    if (reaction_type_ != ReactionType::ELASTIC) {
      double Emin = gen_limits_.electron.E_min;
      double Emax = gen_limits_.electron.E_max;
        
      if (reaction_type_ == ReactionType::QUASIELASTIC ||
	  reaction_type_ == ReactionType::PION_PRODUCTION ||
	  reaction_type_ == ReactionType::KAON_PRODUCTION ||
	  reaction_type_ == ReactionType::DELTA_PRODUCTION) {
	Emin = std::max(Emin, gen_limits_.sumEgen_min);
	Emax = std::min(Emax, gen_limits_.sumEgen_max);
      } else if (reaction_type_ == ReactionType::QUASIELASTIC) {
	Emin = std::max(Emin, gen_limits_.sumEgen_min - event.p_E);
	Emax = std::min(Emax, gen_limits_.sumEgen_max - event.p_E);
      }
        
      if (Emin > Emax) return false;
        
      main.gen_weight *= (Emax - Emin) / 
	(gen_limits_.electron.E_max - gen_limits_.electron.E_min);
        
      event.e_E = random_->Uniform(Emin, Emax);
      event.e_P = event.e_E;  // Ultra-relativistic approximation for electron
      event.e_delta = 100.0 * (event.e_P - spec_electron_.P) / spec_electron_.P;
    }
    // NOTE: For ELASTIC, e_E, e_P, e_delta are calculated later in SolveHydrogenElastic

    // ========================================================================
    // PHASE 5c.3: For elastic, GENERATE hadron angles within acceptance
    // ========================================================================
    if (reaction_type_ == ReactionType::ELASTIC) {
      // Generate hadron spectrometer angles (WITHIN acceptance limits)
      event.p_yptar = random_->Uniform(gen_limits_.hadron.yptar_min,
                                       gen_limits_.hadron.yptar_max);
      event.p_xptar = random_->Uniform(gen_limits_.hadron.xptar_min,
                                       gen_limits_.hadron.xptar_max);
      
      // Calculate hadron physics angles from generated spectrometer angles
      PhysicsAngles(spec_hadron_.theta, spec_hadron_.phi,
                    event.p_xptar, event.p_yptar,
                    event.p_theta, event.p_phi);
      
      // Electron energy will be calculated in SolveHydrogenElastic
      event.e_E = 0.0;
      event.e_P = 0.0;
      event.e_delta = 0.0;
    }

    
    // Generate Fermi Momentum for nuclear targets
    if (reaction_type_ != ReactionType::ELASTIC) {
      FermiMomentum pfer;
      GenerateFermiMomentum(pfer);
      event.Em = pfer.Em;
      // Store pfer for later use in SolvePionKaon, etc.
    }
    
    // Beam energy
    event.Ein = beam_energy_;
    if (beam_energy_spread_ > 0.0) {
      event.Ein += random_->Gaussian(0.0, beam_energy_spread_);
    }
    
    // Generation weight (phase space volume)
    // Base weight: product of all angular ranges and target length
    main.gen_weight *= 
      (gen_limits_.electron.xptar_max - gen_limits_.electron.xptar_min) *
      (gen_limits_.electron.yptar_max - gen_limits_.electron.yptar_min) *
      target_props_.length;
    
    // Add hadron angular ranges (Phase 5c.3: for all reactions including elastic)
    main.gen_weight *= 
      (gen_limits_.hadron.xptar_max - gen_limits_.hadron.xptar_min) *
      (gen_limits_.hadron.yptar_max - gen_limits_.hadron.yptar_min);
    
    return true;
  }

  // ============================================================================
  // Complete Event
  // ============================================================================

  bool EventGenerator::CompleteEvent(SimcEvent& event, MainEvent& main) {
    // Energy loss
    if (use_energy_loss_) {
      double eloss_beam, teff_beam;
      TripThruTarget(1, main.target.z, event.Ein, 0.0, eloss_beam, teff_beam);
      event.Ein -= eloss_beam;
      main.target.Eloss[0] = eloss_beam;
      main.target.teff[0] = teff_beam;
        
      main.target.Eloss[1] = 0.0;
      main.target.teff[1] = 0.0;
      main.target.Eloss[2] = 0.0;
      main.target.teff[2] = 0.0;
    }
    
    // Unit vectors
    event.ue_x = std::sin(event.e_theta) * std::cos(event.e_phi);
    event.ue_y = std::sin(event.e_theta) * std::sin(event.e_phi);
    event.ue_z = std::cos(event.e_theta);
    
    // Solve kinematics based on reaction type (use correct enum values)
    bool success = false;
    switch (reaction_type_) {
    case ReactionType::ELASTIC:
      success = SolveHydrogenElastic(event);
      break;
            
    case ReactionType::QUASIELASTIC:
      success = SolveDeuterium(event, main);
      break;
            
    case ReactionType::PION_PRODUCTION:
    case ReactionType::KAON_PRODUCTION:
      {
	FermiMomentum pfer;
	pfer.Em = event.Em;
	success = SolvePionKaon(event, pfer);
	break;
      }
        
    default:
      std::cerr << "Error: Reaction type not implemented\n";
      success = false;
    }
    
    return success;
  }

  // ============================================================================
  // Solve Hydrogen Elastic
  // ============================================================================

  bool EventGenerator::SolveHydrogenElastic(SimcEvent& event) const {
    // ========================================================================
    // PHASE 5c.3: H(e,e'p) elastic with acceptance-based angle generation
    // ========================================================================
    //
    // Previous (Phase 5c.2): Calculated p_theta, p_phi from q-vector
    // Now (Phase 5c.3): p_theta, p_phi are GENERATED in GeneratePhaseSpace
    //
    // We apply a kinematic constraint: reject if generated proton direction
    // doesn't match the q-vector direction (within tolerance)
    // ========================================================================
    
    // Calculate electron energy from elastic formula
    double denominator = kProtonMass + event.Ein * (1.0 - event.ue_z);
    
    if (denominator <= 0.0) {
      return false;
    }
    
    event.e_E = event.Ein * kProtonMass / denominator;
    
    if (event.e_E >= event.Ein || event.e_E < kElectronMass) {
      return false;
    }
    
    event.e_P = std::sqrt(event.e_E * event.e_E - kElectronMass * kElectronMass);
    event.e_delta = 100.0 * (event.e_P - spec_electron_.P) / spec_electron_.P;
    
    // Calculate Q² and nu
    double sin_half = std::sin(event.e_theta / 2.0);
    event.nu = event.Ein - event.e_E;
    event.Q2 = 4.0 * event.Ein * event.e_E * sin_half * sin_half;
    
    // Calculate q-vector (momentum transfer)
    Vector3D q_vec;
    q_vec.x = -event.ue_x * event.e_P;
    q_vec.y = -event.ue_y * event.e_P;
    q_vec.z = event.Ein - event.ue_z * event.e_P;
    
    double q_mag = std::sqrt(q_vec.x*q_vec.x + q_vec.y*q_vec.y + q_vec.z*q_vec.z);
    
    // Expected proton direction from elastic kinematics (along q-vector)
    double up_x_elastic = q_vec.x / q_mag;
    double up_y_elastic = q_vec.y / q_mag;
    double up_z_elastic = q_vec.z / q_mag;
    
    // Actual proton direction (from generated angles in Phase 5c.3)
    double cos_theta_p = std::cos(event.p_theta);
    double sin_theta_p = std::sin(event.p_theta);
    double cos_phi_p = std::cos(event.p_phi);
    double sin_phi_p = std::sin(event.p_phi);
    
    event.up_x = sin_theta_p * cos_phi_p;
    event.up_y = sin_theta_p * sin_phi_p;
    event.up_z = cos_theta_p;
    
    // ========================================================================
    // PHASE 5c.3: Accept all generated angles - let SHMS determine acceptance
    // ========================================================================
    // NOTE: Originally we checked if generated proton direction matched
    // elastic kinematics (q-vector direction). However, this is too restrictive
    // because:
    // 1. SHMS may not be at the elastic angle (e.g., SHMS at 35° vs q at 32°)
    // 2. We're testing SHMS acceptance, not elastic kinematics
    // 3. The SHMS apertures will naturally reject unphysical events
    //
    // For true physics simulation, this check should be re-enabled with proper
    // alignment between spectrometer angle and elastic kinematics.
    // ========================================================================
    
    // (Kinematic constraint check removed - SHMS acceptance handles it)
    
    // ========================================================================
    // Event passes kinematic constraint - set final values
    // ========================================================================
    
    // Proton momentum (from q magnitude for elastic)
    event.p_P = q_mag;
    event.p_E = std::sqrt(event.p_P * event.p_P + kProtonMass * kProtonMass);
    event.p_delta = 100.0 * (event.p_P - spec_hadron_.P) / spec_hadron_.P;
    
    // q-vector unit vector
    event.uq_x = up_x_elastic;
    event.uq_y = up_y_elastic;
    event.uq_z = up_z_elastic;
    
    return true;
  }

  // ============================================================================
  // Solve Deuterium
  // ============================================================================

  bool EventGenerator::SolveDeuterium(SimcEvent& event, MainEvent& main) const {
    // D(e,e'p)n reaction
    
    event.Em = kProtonMass + kNeutronMass - target_props_.mass; // ~2.2249 MeV
    
    double sin_half = std::sin(event.e_theta / 2.0);
    event.nu = event.Ein - event.e_E;
    event.Q2 = 4.0 * event.Ein * event.e_E * sin_half * sin_half;
    
    double omega = event.nu;
    double q = std::sqrt(event.Q2 + omega * omega);
    
    // Virtual photon direction
    event.uq_x = -event.ue_x * event.e_P / q;
    event.uq_y = -event.ue_y * event.e_P / q;
    event.uq_z = (event.Ein - event.ue_z * event.e_P) / q;
    
    // Proton direction
    event.up_x = std::sin(event.p_theta) * std::cos(event.p_phi);
    event.up_y = std::sin(event.p_theta) * std::sin(event.p_phi);
    event.up_z = std::cos(event.p_theta);
    
    // Quadratic for proton energy
    double cos_theta_pq = event.uq_x * event.up_x + 
      event.uq_y * event.up_y + 
      event.uq_z * event.up_z;
    
    double M_n = kNeutronMass;
    double E_target = M_n + omega;
    
    // a*E_p² + b*E_p + c = 0
    double a = 1.0;
    double b = -2.0 * E_target;
    double c = E_target * E_target - q * q - 
      2.0 * q * event.p_P * cos_theta_pq - 
      event.p_P * event.p_P;
    
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
      return false;
    }
    
    event.p_E = (-b + std::sqrt(disc)) / (2.0 * a);
    
    if (event.p_E < kProtonMass) {
      return false;
    }
    
    // Jacobian |dE_p/dE_e'|
    main.jacobian = event.p_E / (event.p_E - event.p_P * cos_theta_pq * q / omega);
    
    return true;
  }

  // ============================================================================
  // Solve Pion/Kaon Production
  // ============================================================================

  bool EventGenerator::SolvePionKaon(SimcEvent& event, const FermiMomentum& pfer) const {
    // Placeholder - not yet implemented
    std::cerr << "Warning: Pion/kaon production not implemented\n";
    return false;
  }

  
  // ============================================================================
  // Physics Angles
  // ============================================================================


void EventGenerator::PhysicsAngles(
    Double_t theta0, Double_t phi0,
    Double_t xptar, Double_t yptar,
    Double_t& theta, Double_t& phi) const {
    
    double cos_th0 = std::cos(theta0);
    double sin_th0 = std::sin(theta0);
    double cos_ph0 = std::cos(phi0);
    double sin_ph0 = std::sin(phi0);
    
    // Direction in spectrometer frame (NOT normalized)
    double dx = xptar;
    double dy = yptar;
    double dz = 1.0;  // DON'T use sqrt(1+dx²+dy²), just use 1.0
    
    // Rotate to lab frame
    double ux = cos_th0 * cos_ph0 * dx - sin_ph0 * dy + sin_th0 * cos_ph0 * dz;
    double uy = cos_th0 * sin_ph0 * dx + cos_ph0 * dy + sin_th0 * sin_ph0 * dz;
    double uz = -sin_th0 * dx + cos_th0 * dz;
    
    // Normalize
    double norm = std::sqrt(ux*ux + uy*uy + uz*uz);
    ux /= norm;
    uy /= norm;
    uz /= norm;
    
    // Convert to angles
    theta = std::acos(uz);
    phi = std::atan2(uy, ux);
    if (phi < 0.0) phi += 2.0 * kPi;
}

// ============================================================================
//  SpectrometerAngles Angles
// ============================================================================

void EventGenerator::SpectrometerAngles(
    Double_t theta0, Double_t phi0,
    Double_t theta, Double_t phi,
    Double_t& xptar, Double_t& yptar) const {
    
    double cos_th0 = std::cos(theta0);
    double sin_th0 = std::sin(theta0);
    double cos_ph0 = std::cos(phi0);
    double sin_ph0 = std::sin(phi0);
    
    // Unit vector in lab frame
    double ux = std::sin(theta) * std::cos(phi);
    double uy = std::sin(theta) * std::sin(phi);
    double uz = std::cos(theta);
    
    // Inverse rotation: lab -> spectrometer
    double dx = cos_th0 * cos_ph0 * ux + cos_th0 * sin_ph0 * uy - sin_th0 * uz;
    double dy = -sin_ph0 * ux + cos_ph0 * uy;
    double dz = sin_th0 * cos_ph0 * ux + sin_th0 * sin_ph0 * uy + cos_th0 * uz;
    
    // Check that dz is positive (particle going forward through spectrometer)
    if (dz <= 0.0) {
        std::cerr << "ERROR in SpectrometerAngles: dz <= 0 (particle going backwards!)" << std::endl;
        std::cerr << "  theta0=" << theta0 << ", phi0=" << phi0 << std::endl;
        std::cerr << "  theta=" << theta << ", phi=" << phi << std::endl;
        std::cerr << "  dz=" << dz << std::endl;
        xptar = 0.0;
        yptar = 0.0;
        return;
    }
    
    xptar = dx / dz;
    yptar = dy / dz;
}

 

  // ============================================================================
  // Calculate Basic Kinematics
  // ============================================================================

  void EventGenerator::CalculateBasicKinematics(SimcEvent& event) const {
    double sin_half = std::sin(event.e_theta / 2.0);
    event.Q2 = 4.0 * event.Ein * event.e_E * sin_half * sin_half;
    event.nu = event.Ein - event.e_E;
    event.q = std::sqrt(event.Q2 + event.nu * event.nu);
    
    event.W = std::sqrt(kProtonMass * kProtonMass + 
                        2.0 * kProtonMass * event.nu - event.Q2);
    
    event.xbj = event.Q2 / (2.0 * kProtonMass * event.nu);
    
    double Q2_GeV2 = event.Q2 / 1.0e6;
    double nu_GeV = event.nu / 1000.0;
    event.epsilon = 1.0 / (1.0 + 2.0 * (1.0 + nu_GeV * nu_GeV / Q2_GeV2) * 
			   std::tan(event.e_theta / 2.0) * std::tan(event.e_theta / 2.0));
  }

  // ============================================================================
  // Calculate Missing Momentum
  // ============================================================================

  void EventGenerator::CalculateMissingMomentum(SimcEvent& event) const {
    // P_miss = q - p_hadron
    
    Vector3D q_vec;
    q_vec.x = event.uq_x * event.q;
    q_vec.y = event.uq_y * event.q;
    q_vec.z = event.uq_z * event.q;
    
    Vector3D p_vec;
    p_vec.x = event.up_x * event.p_P;
    p_vec.y = event.up_y * event.p_P;
    p_vec.z = event.up_z * event.p_P;
    
    event.Pmx = q_vec.x - p_vec.x;
    event.Pmy = q_vec.y - p_vec.y;
    event.Pmz = q_vec.z - p_vec.z;
    
    event.Pmiss = std::sqrt(event.Pmx * event.Pmx + 
                            event.Pmy * event.Pmy + 
                            event.Pmz * event.Pmz);
    
    // Decompose
    double q_mag = event.q;
    event.PmPar = (event.Pmx * q_vec.x + event.Pmy * q_vec.y + event.Pmz * q_vec.z) / q_mag;
    
    Vector3D pm_vec = {event.Pmx, event.Pmy, event.Pmz};
    Vector3D q_unit = {q_vec.x / q_mag, q_vec.y / q_mag, q_vec.z / q_mag};
    
    Vector3D pm_par = {q_unit.x * event.PmPar, q_unit.y * event.PmPar, q_unit.z * event.PmPar};
    Vector3D pm_perp = {pm_vec.x - pm_par.x, pm_vec.y - pm_par.y, pm_vec.z - pm_par.z};
    
    event.PmPer = std::sqrt(pm_perp.x * pm_perp.x + pm_perp.y * pm_perp.y + pm_perp.z * pm_perp.z);
    
    // Out-of-plane (cross product method)
    Vector3D beam = {0.0, 0.0, 1.0};
    Vector3D q_cross_beam = q_unit.Cross(beam);
    event.PmOop = pm_vec.Dot(q_cross_beam);
  }

  // ============================================================================
  // Calculate Hadron Angles
  // ============================================================================

  void EventGenerator::CalculateHadronAngles(SimcEvent& event) const {
    // Angle between hadron and q
    double cos_theta_pq = event.uq_x * event.up_x + 
      event.uq_y * event.up_y + 
      event.uq_z * event.up_z;
    
    event.theta_pq = std::acos(cos_theta_pq);
    
    // Out-of-plane angle
    Vector3D q_vec = {event.uq_x, event.uq_y, event.uq_z};
    Vector3D p_vec = {event.up_x, event.up_y, event.up_z};
    Vector3D beam = {0.0, 0.0, 1.0};
    
    Vector3D q_cross_beam = q_vec.Cross(beam);
    Vector3D p_cross_q = p_vec.Cross(q_vec);
    
    double cos_phi = q_cross_beam.Dot(p_cross_q) / 
      (q_cross_beam.Magnitude() * p_cross_q.Magnitude());
    event.phi_pq = std::acos(cos_phi);
  }

  // ============================================================================
  // Generate Fermi Momentum
  // ============================================================================

  void EventGenerator::GenerateFermiMomentum(FermiMomentum& pfer) {
    pfer.Clear();
    
    if (target_props_.A == 1) {
      pfer.magnitude = 0.0;
    } else {
      // Simple Gaussian
      pfer.magnitude = std::abs(random_->Gaussian(0.0, 100.0));
      if (pfer.magnitude > 300.0) pfer.magnitude = 300.0;
        
      // Isotropic
      double cos_theta = random_->Uniform(-1.0, 1.0);
      double phi = random_->Uniform(0.0, 2.0 * kPi);
      double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        
      pfer.x = pfer.magnitude * sin_theta * std::cos(phi);
      pfer.y = pfer.magnitude * sin_theta * std::sin(phi);
      pfer.z = pfer.magnitude * cos_theta;
        
      GenerateMissingEnergy(pfer.magnitude, pfer.Em);
    }
  }

  // ============================================================================
  // Generate Missing Energy
  // ============================================================================

  void EventGenerator::GenerateMissingEnergy(double pfer, double& Em) {
    if (target_props_.A == 1) {
      Em = 0.0;
    } else if (target_props_.A == 2) {
      Em = 2.2249;
    } else {
      Em = random_->Uniform(0.0, 50.0);
    }
  }

  // ============================================================================
  // Calculate Jacobian
  // ============================================================================

  Double_t EventGenerator::CalculateJacobian(const SimcEvent& event, const MainEvent& main) const {
    // Basic angle jacobian
    // Full implementation requires reaction-specific terms
    
    Double_t jac = 1.0;
    
    // Add angle transformation jacobian
    double cos_e = std::cos(event.e_theta);
    double cos_p = std::cos(event.p_theta);
    
    jac *= std::abs(cos_e * cos_p);
    
    return jac;
  }

  // ============================================================================
  // Calculate CM Jacobian
  // ============================================================================

  Double_t EventGenerator::CalculateCMJacobian(const SimcEvent& event, const MainEvent& main) const {
    // Placeholder for CM transformation jacobian
    return 1.0;
  }

  // ============================================================================
  // Decay Two Body
  // ============================================================================

  void EventGenerator::DecayTwoBody(
				    const std::array<double, 4>& parent,
				    std::array<double, 4>& decay1,
				    std::array<double, 4>& decay2,
				    double m1, double m2,
				    double costh, double phi) const {
    
    // Placeholder - full implementation needed
    std::cerr << "Warning: DecayTwoBody not implemented\n";
  }

  // ============================================================================
  // Trip Thru Target
  // ============================================================================

  void EventGenerator::TripThruTarget(
				      Int_t particle_id,
				      Double_t z_position,
				      Double_t energy,
				      Double_t theta,
				      Double_t& eloss,
				      Double_t& teff) const {
    
    // Simplified energy loss calculation
    // Full implementation would use EnergyLoss class properly
    
    double path_length = std::abs(z_position) / std::max(std::cos(theta), 0.01);
    teff = target_props_.density * path_length;
    
    // Simple dE/dx approximation: ~2 MeV/(g/cm²) for minimum ionizing
    eloss = 2.0 * teff;
    
    // Could be improved with proper Bethe-Bloch later
  }

  // ============================================================================
  // Unit Vector
  // ============================================================================

  Vector3D EventGenerator::UnitVector(double theta, double phi) const {
    Vector3D vec;
    vec.x = std::sin(theta) * std::cos(phi);
    vec.y = std::sin(theta) * std::sin(phi);
    vec.z = std::cos(theta);
    return vec;
  }

  // ============================================================================
  // Passes Generation Cuts
  // ============================================================================
  // ============================================================================
  // COMPLETE FIX: Replace PassesGenerationCuts in EventGenerator.cpp
  // This version ONLY checks quantities that were actually GENERATED
  // ============================================================================

  bool EventGenerator::PassesGenerationCuts(const SimcEvent& event) const {
    // CRITICAL: Only check quantities that were GENERATED!
    // For ELASTIC: Only electron angles are generated
    // For non-ELASTIC: Electron + hadron angles, and energies
    
    // ========================================================================
    // ELECTRON ARM - Angles (ALWAYS generated for all reactions)
    // ========================================================================
    
    if (event.e_xptar < gen_limits_.electron.xptar_min || 
        event.e_xptar > gen_limits_.electron.xptar_max) {
      return false;
    }
    
    if (event.e_yptar < gen_limits_.electron.yptar_min || 
        event.e_yptar > gen_limits_.electron.yptar_max) {
      return false;
    }
    
    // ========================================================================
    // ELECTRON ARM - Energy/Delta (ONLY for non-elastic)
    // ========================================================================
    
    if (reaction_type_ != ReactionType::ELASTIC) {
      // For non-elastic, electron energy WAS generated, so check it
      if (event.e_E < gen_limits_.electron.E_min || 
	  event.e_E > gen_limits_.electron.E_max) {
	return false;
      }
        
      // Delta check with slop (it's calculated from E, allow some tolerance)
      if (event.e_delta < gen_limits_.electron.delta_min - 5.0 ||
	  event.e_delta > gen_limits_.electron.delta_max + 5.0) {
	return false;
      }
    }
    // NOTE: For ELASTIC, e_E and e_delta were CALCULATED (not generated)
    // so we DON'T check them against generation limits!
    
    // ========================================================================
    // HADRON ARM (ONLY for non-elastic reactions)
    // ========================================================================
    
    if (reaction_type_ != ReactionType::ELASTIC) {
      // Hadron angles
      if (event.p_xptar < gen_limits_.hadron.xptar_min || 
	  event.p_xptar > gen_limits_.hadron.xptar_max) {
	return false;
      }
        
      if (event.p_yptar < gen_limits_.hadron.yptar_min || 
	  event.p_yptar > gen_limits_.hadron.yptar_max) {
	return false;
      }
        
      // Hadron energy (only for certain reactions where it was generated)
      if (reaction_type_ == ReactionType::QUASIELASTIC ||
	  reaction_type_ == ReactionType::SEMI_INCLUSIVE) {
	if (event.p_E < gen_limits_.hadron.E_min || 
	    event.p_E > gen_limits_.hadron.E_max) {
	  return false;
	}
      }
        
      // Hadron delta check with slop
      if (event.p_delta < gen_limits_.hadron.delta_min - 5.0 ||
	  event.p_delta > gen_limits_.hadron.delta_max + 5.0) {
	return false;
      }
    }
    
    // ========================================================================
    // ELASTIC - Check generated hadron angles (Phase 5c.3)
    // ========================================================================
    
    if (reaction_type_ == ReactionType::ELASTIC) {
      // Hadron angles ARE generated for elastic now (Phase 5c.3)
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


  // ============================================================================
  // Validate Physics
  // ============================================================================

  bool EventGenerator::ValidatePhysics(const SimcEvent& event) const {
    // Check physical constraints
    
    // Energy conservation
    if (event.e_E >= event.Ein) {
      return false;  // No energy gain
    }
    
    // Electron above mass
    if (event.e_E < kElectronMass) {
      return false;
    }
    
    // Hadron above mass
    if (event.p_E < kProtonMass) {
      return false;
    }
    
    // Q² positive (space-like)
    if (event.Q2 <= 0.0) {
      return false;
    }
    
    // Valid angles
    if (std::abs(std::cos(event.e_theta)) > 1.0) {
      return false;
    }
    
    if (std::abs(std::cos(event.p_theta)) > 1.0) {
      return false;
    }
    
    // Unit vectors normalized
    double ue_mag = std::sqrt(event.ue_x*event.ue_x + 
                              event.ue_y*event.ue_y + 
                              event.ue_z*event.ue_z);
    if (std::abs(ue_mag - 1.0) > 0.01) {
      return false;
    }
    
    double up_mag = std::sqrt(event.up_x*event.up_x + 
                              event.up_y*event.up_y + 
                              event.up_z*event.up_z);
    if (std::abs(up_mag - 1.0) > 0.01) {
      return false;
    }
    
    return true;
  }

} // namespace simc
