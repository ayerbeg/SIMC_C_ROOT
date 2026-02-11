// src/physics/EnergyLoss.cpp
// Implementation of energy loss calculations
// Ported from enerloss_new.f

#include "simc/physics/EnergyLoss.h"
#include "simc/core/SimcConstants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cctype>

namespace simc {
  using namespace constants;

  // Static error counter for warnings
  static int eloss_error_count = 0;

  // ============================================================================
  // Main Energy Loss Calculation
  // ============================================================================

  double EnergyLoss::Calculate(double length, double density, 
			       double Z, double A,
			       double E_particle, double M_particle,
			       RandomGenerator& rng,
			       LossType type) {
    // Port from enerloss_new() in enerloss_new.f
    
    // Calculate thickness in g/cm²
    double thickness = length * density;
    
    if (thickness <= 0.0) {
      return 0.0;
    }
    
    // Check for valid particle energy
    if (E_particle <= M_particle) {
      std::cerr << "Warning: E_particle <= M_particle, setting Eloss = 0\n";
      return 0.0;
    }
    
    // Calculate gamma and beta
    double gamma = E_particle / M_particle;
    
    // Safety check for gamma
    if (gamma < 1.0) {
      gamma = 1.0 + 1e-10;  // Ensure gamma >= 1
    }
    
    double beta2 = 1.0 - 1.0 / (gamma * gamma);
    
    // Safety check for beta
    if (beta2 <= 0.0 || beta2 >= 1.0) {
      std::cerr << "Warning: Invalid beta^2 = " << beta2 << ", using beta = 0.9\n";
      beta2 = 0.81;  // beta = 0.9
    }
    
    double beta = std::sqrt(beta2);
    
    // Get ionization potential
    double I = GetIonizationPotential(Z);
    
    // Calculate most probable energy loss
    double Eloss_mp = MostProbableEnergyLoss(thickness, Z, A, beta, gamma, density);
    
    // Calculate xi parameter for Landau fluctuations
    double xi = CalcXi(thickness, Z, A, beta);
    
    // Calculate lambda based on loss type
    double lambda_val;
    switch (type) {
    case LossType::SAMPLED: {
      // Sample from Landau distribution
      lambda_val = SampleLandau(rng);
      break;
    }
    case LossType::MINIMUM:
      lambda_val = 3.0;
      break;
    case LossType::MAXIMUM:
      lambda_val = 0.0067;
      break;
    case LossType::MOST_PROBABLE:
      lambda_val = 1.0;
      break;
    default:
      lambda_val = 1.0;
    }
    
    // Total energy loss = Landau fluctuation + most probable
    double Eloss = lambda_val * xi + Eloss_mp;
    
    // Safety check for negative energy loss
    if (Eloss < 0.0) {
      Eloss = 0.0;
    }
    
    // Check that energy loss doesn't exceed kinetic energy
    double KE = E_particle - M_particle;
    if (Eloss > KE) {
      Eloss = KE - 0.0000001;  // Leave tiny amount
        
      ++eloss_error_count;
      if (eloss_error_count <= 10) {
	std::cerr << "Warning: Eloss > Total KE; forcing Eloss = KE" << std::endl;
	if (eloss_error_count == 10) {
	  std::cerr << "  FURTHER ELOSS ERRORS SUPPRESSED" << std::endl;
	}
      }
    }
    
    return Eloss;
  }

  // ============================================================================
  // Most Probable Energy Loss (Bethe-Bloch with corrections)
  // ============================================================================

  double EnergyLoss::MostProbableEnergyLoss(double thickness, 
					    double Z, double A,
					    double beta, double gamma,
					    double density) {
    // Port from enerloss_new.f
    // Modern Bethe-Bloch formula
    
    if (thickness <= 0.0 || beta <= 0.0 || beta >= 1.0) {
      return 0.0;
    }
    
    double I = GetIonizationPotential(Z);
    double denscorr = DensityCorrection(beta, gamma, Z, A, density);
    
    double beta2 = beta * beta;
    
    // Safety check for log arguments
    double log_arg1 = Me / (I * I);
    if (log_arg1 <= 0.0) {
      log_arg1 = 1.0;
    }
    
    double log_arg2 = 0.1536 * Z / A * thickness / beta2;
    if (log_arg2 <= 0.0) {
      log_arg2 = 1.0;
    }
    
    // Bethe-Bloch formula (in MeV)
    // 0.1536 MeV cm²/g is the constant for dE/dx
    double Eloss_mp = 0.1536e-03 * Z / A * thickness / beta2 * (
								std::log(log_arg1) + 1.063 + 2.0 * std::log(gamma * beta) +
								std::log(log_arg2) - beta2 - denscorr
								);
    
    // Convert to MeV (formula gives GeV, we multiply by 1000)
    double result = Eloss_mp * 1000.0;
    
    // Safety check for unreasonable values
    if (result < 0.0) {
      result = 0.0;
    }
    if (result > 1000.0) {  // More than 1 GeV is suspicious
      std::cerr << "Warning: Very large energy loss calculated: " << result << " MeV\n";
    }
    
    return result;
  }

  // ============================================================================
  // Bethe-Bloch dE/dx
  // ============================================================================

  double EnergyLoss::BetheBloch(double beta, double gamma,
				double Z, double A, double I,
				double density) {
    // Bethe-Bloch formula for dE/dx
    // Returns MeV/(g/cm²)
    
    if (beta <= 0.0 || beta >= 1.0) {
      return 0.0;
    }
    
    if (gamma < 1.0) {
      gamma = 1.0;
    }
    
    double Tmax = 2.0 * Me * beta * beta * gamma * gamma;
    double denscorr = DensityCorrection(beta, gamma, Z, A, density);
    
    // Bethe-Bloch constant (MeV cm²/g)
    const double K = 0.307075;  // MeV cm²/g
    
    double beta2 = beta * beta;
    
    // Safety check for log argument
    double log_arg = 2.0 * Me * beta2 * gamma * gamma * Tmax / (I * I);
    if (log_arg <= 0.0) {
      log_arg = 1.0;
    }
    
    double dEdx = K * Z / A / beta2 * (
				       0.5 * std::log(log_arg) - beta2 - denscorr / 2.0
				       );
    
    return dEdx;
  }

  // ============================================================================
  // Ionization Potential
  // ============================================================================

  double EnergyLoss::GetIonizationPotential(double Z) {
    // Port from enerloss_new.f
    
    if (Z <= 0.0) {
      return 21.8e-06;  // Default to hydrogen
    }
    
    if (Z == 1.0) {
      // Hydrogen
      return 21.8e-06;  // MeV (21.8 eV)
    } else {
      // Other elements: I = 16*Z^0.9 eV
      return 16.0 * std::pow(Z, 0.9) * 1.0e-06;  // MeV
    }
  }

  // ============================================================================
  // Density Effect Correction
  // ============================================================================

  double EnergyLoss::DensityCorrection(double beta, double gamma,
				       double Z, double A, double density) {
    // Port from enerloss_new.f
    // Density effect correction (from J. Volmer)
    
    if (beta <= 0.0 || gamma < 1.0) {
      return 0.0;
    }
    
    double I = GetIonizationPotential(Z);
    double hnup = PlasmaFrequency(density, Z, A);
    
    // Safety checks for log arguments
    if (hnup <= 0.0) hnup = 1e-10;
    if (I <= 0.0) I = 1e-10;
    
    double log10bg = std::log10(beta * gamma);
    double C0 = std::log(hnup) - std::log(I) + 0.5;
    
    double denscorr;
    
    if (log10bg < 0.0) {
      denscorr = 0.0;
    } else if (log10bg < 3.0) {
      denscorr = C0 + std::log(10.0) * log10bg + 
	std::abs(C0 / 27.0) * std::pow(3.0 - log10bg, 3);
    } else if (log10bg < 4.7) {
      denscorr = C0 + std::log(10.0) * log10bg;
    } else {
      denscorr = C0 + std::log(10.0) * 4.7;
    }
    
    return denscorr;
  }

  // ============================================================================
  // Plasma Frequency
  // ============================================================================

  double EnergyLoss::PlasmaFrequency(double density, double Z, double A) {
    // Port from enerloss_new.f
    // Returns hbar*omega_p in MeV
    
    if (density <= 0.0 || Z <= 0.0 || A <= 0.0) {
      return 1e-10;  // Small positive value to avoid log(0)
    }
    
    return 28.816e-06 * std::sqrt(density * Z / A);
  }

  // ============================================================================
  // Xi Parameter (Landau)
  // ============================================================================

  double EnergyLoss::CalcXi(double thickness, double Z, double A, double beta) {
    // Port from enerloss_new.f
    // Xi parameter for Landau distribution
    
    if (thickness <= 0.0 || beta <= 0.0 || beta >= 1.0) {
      return 0.0;
    }
    
    double beta2 = beta * beta;
    return 0.307075 / 2.0 * Z / A * thickness / beta2;
  }

  // ============================================================================
  // Sample from Landau Distribution
  // ============================================================================

  double EnergyLoss::SampleLandau(RandomGenerator& rng) {
    // Sample lambda = -2*ln(x) where x is from Gaussian tail
    // Use exponential tail approximation for large nsigma
    
    double x = GaussianTail(rng, 10.0);
    
    if (x > 0.0) {
      return -2.0 * std::log(x);
    } else {
      return 100000.0;  // Very large value
    }
  }

  // ============================================================================
  // Gaussian Tail Sampling
  // ============================================================================

  double EnergyLoss::GaussianTail(RandomGenerator& rng, double nsigma) {
    // Sample from Gaussian tail beyond nsigma standard deviations
    // For large nsigma (>3), use exponential approximation
    // For small nsigma, use rejection sampling
    
    if (nsigma > 3.0) {
      // Exponential approximation for extreme tail
      // Gaussian tail: f(x) ~ exp(-x²/2) for x >> nsigma
      // Approximate as: f(x) ~ exp(-nsigma*(x-nsigma)) for x > nsigma
      double rate = nsigma;
      double tail = rng.Exponential(rate);
      return nsigma + tail;
    } else {
      // Rejection sampling for moderate nsigma
      const int MAX_TRIES = 10000;
      int tries = 0;
        
      while (tries < MAX_TRIES) {
	double x = std::abs(rng.Gaussian(0.0, 1.0));
	if (x >= nsigma) {
	  return x;
	}
	++tries;
      }
        
      // Fallback after max tries (should be extremely rare)
      std::cerr << "Warning: GaussianTail exceeded max tries, using exponential fallback\n";
      return nsigma + rng.Exponential(1.0);
    }
  }

  // ============================================================================
  // Material Database
  // ============================================================================

  Material Material::LiquidHydrogen() {
    Material mat;
    mat.name = "LH2";
    mat.Z = 1.0;
    mat.A = 1.008;
    mat.density = 0.0708;  // g/cm³ at 20K
    mat.X0 = 890.0;        // g/cm²
    return mat;
  }

  Material Material::LiquidDeuterium() {
    Material mat;
    mat.name = "LD2";
    mat.Z = 1.0;
    mat.A = 2.014;
    mat.density = 0.162;   // g/cm³
    mat.X0 = 1260.0;       // g/cm²
    return mat;
  }

  Material Material::Helium3() {
    Material mat;
    mat.name = "He3";
    mat.Z = 2.0;
    mat.A = 3.016;
    mat.density = 0.059;   // g/cm³ (liquid)
    mat.X0 = 755.0;        // g/cm²
    return mat;
  }

  Material Material::Helium4() {
    Material mat;
    mat.name = "He4";
    mat.Z = 2.0;
    mat.A = 4.003;
    mat.density = 0.125;   // g/cm³ (liquid)
    mat.X0 = 755.0;        // g/cm²
    return mat;
  }

  Material Material::Carbon() {
    Material mat;
    mat.name = "C";
    mat.Z = 6.0;
    mat.A = 12.011;
    mat.density = 2.26;    // g/cm³ (graphite)
    mat.X0 = 42.7;         // g/cm²
    return mat;
  }

  Material Material::Aluminum() {
    Material mat;
    mat.name = "Al";
    mat.Z = 13.0;
    mat.A = 26.982;
    mat.density = 2.70;    // g/cm³
    mat.X0 = 24.01;        // g/cm²
    return mat;
  }

  Material Material::Iron() {
    Material mat;
    mat.name = "Fe";
    mat.Z = 26.0;
    mat.A = 55.845;
    mat.density = 7.87;    // g/cm³
    mat.X0 = 13.84;        // g/cm²
    return mat;
  }

  Material Material::Copper() {
    Material mat;
    mat.name = "Cu";
    mat.Z = 29.0;
    mat.A = 63.546;
    mat.density = 8.96;    // g/cm³
    mat.X0 = 12.86;        // g/cm²
    return mat;
  }

  Material Material::Gold() {
    Material mat;
    mat.name = "Au";
    mat.Z = 79.0;
    mat.A = 196.967;
    mat.density = 19.32;   // g/cm³
    mat.X0 = 6.46;         // g/cm²
    return mat;
  }

  Material Material::GetByName(const std::string& name) {
    // Create lowercase copy for case-insensitive comparison
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), 
                   [](unsigned char c){ return std::tolower(c); });
    
    if (lower_name == "lh2" || lower_name == "h2") return LiquidHydrogen();
    if (lower_name == "ld2" || lower_name == "d2") return LiquidDeuterium();
    if (lower_name == "he3") return Helium3();
    if (lower_name == "he4") return Helium4();
    if (lower_name == "c" || lower_name == "carbon") return Carbon();
    if (lower_name == "al" || lower_name == "aluminum") return Aluminum();
    if (lower_name == "fe" || lower_name == "iron") return Iron();
    if (lower_name == "cu" || lower_name == "copper") return Copper();
    if (lower_name == "au" || lower_name == "gold") return Gold();
    
    throw std::runtime_error("Unknown material: " + name);
  }

} // namespace simc
