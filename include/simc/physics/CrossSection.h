// include/simc/CrossSection.h
// Cross section calculations for electron scattering
// COMPLETE PORT from physics_proton.f, physics_pion.f, physics_kaon.f

#ifndef SIMC_CROSS_SECTION_H
#define SIMC_CROSS_SECTION_H

#include "simc/core/SimcTypes.h"
#include "simc/core/SimcEvent.h"
#include <memory>

namespace simc {

// Forward declaration
class MainEvent;

/**
 * @class CrossSectionBase
 * @brief Abstract base class for cross section calculations
 * 
 * All reaction-specific cross sections inherit from this class.
 */
class CrossSectionBase {
public:
    virtual ~CrossSectionBase() = default;
    
    /**
     * @brief Calculate differential cross section
     * @param evt Event with kinematics
     * @return Cross section in microbarns/(MeV/c)^n/sr^m
     */
    virtual double Calculate(const SimcEvent& evt) const = 0;
    
    /**
     * @brief Get reaction type
     */
    virtual ReactionType GetReactionType() const = 0;
    
    /**
     * @brief Check if kinematics are physical
     */
    virtual bool IsPhysical(const SimcEvent& evt) const = 0;
};

/**
 * @class ElasticCrossSection
 * @brief Elastic electron-proton scattering
 * 
 * EXACT PORT from physics_proton.f: sigep(), sigMott(), fofa_best_fit()
 * Uses Rosenbluth formula with Peter Bosted form factors
 */
class ElasticCrossSection : public CrossSectionBase {
public:
    ElasticCrossSection() = default;
    
    /**
     * @brief Calculate elastic e-p cross section
     * @param evt Event with kinematics
     * @return Cross section in microbarn/sr
     * 
     * PORT FROM: sigep() in physics_proton.f
     */
    double Calculate(const SimcEvent& evt) const override;
    
    ReactionType GetReactionType() const override { return ReactionType::ELASTIC; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    /**
     * @brief Calculate Mott cross section
     * @param e0 Incident energy (MeV)
     * @param theta Scattering angle (rad)
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @return Mott cross section in microbarn/sr
     * 
     * PORT FROM: sigMott() in physics_proton.f
     */
    static double SigMott(double e0, double theta, double Q2);
    
    /**
     * @brief Get proton form factors using Peter Bosted's fit
     * @param qsquar -Q^2/(hbarc^2) (fm^-2)
     * @param GE Electric form factor (output)
     * @param GM Magnetic form factor (output)
     * 
     * PORT FROM: fofa_best_fit() in physics_proton.f
     * Reference: Phys. Rev. C 51, 409, Eqs. 4 and 5
     */
    static void FofaBestFit(double qsquar, double& GE, double& GM);
};

/**
 * @class QuasiElasticCrossSection
 * @brief Quasi-elastic (e,e'p) scattering from nuclei
 * 
 * EXACT PORT from physics_proton.f: deForest()
 * Includes Fermi motion, binding energy, and off-shell corrections
 */
class QuasiElasticCrossSection : public CrossSectionBase {
public:
    /**
     * @brief Constructor
     * @param deforest_flag Controls which deForest prescription to use:
     *   0 = use sigcc1 (default)
     *   1 = use sigcc2
     *  -1 = use sigcc1 ONSHELL (Ebar = E' - nu, qbar = q)
     */
    explicit QuasiElasticCrossSection(int deforest_flag = 0);
    
    /**
     * @brief Calculate quasi-elastic cross section
     * @param evt Event with kinematics
     * @return Cross section in microbarn*MeV^2/sr^2
     * 
     * PORT FROM: deForest() in physics_proton.f
     * 
     * IMPORTANT NOTES FROM FORTRAN:
     * - Units are microbarn*MeV^2/sr^2 (need to multiply by spectral function S(E,p) in MeV^-4)
     * - Implements de Forest off-shell prescription
     * - Can be called with vertex or recon event records
     * - Formula includes E'/Ebar combined with K=E'*p' to give p'/Ebar
     */
    double Calculate(const SimcEvent& evt) const override;
    
    ReactionType GetReactionType() const override { return ReactionType::QUASIELASTIC; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    /**
     * @brief Set deForest flag
     */
    void SetDeForestFlag(int flag) { deforest_flag_ = flag; }
    
    /**
     * @brief Get deForest flag
     */
    int GetDeForestFlag() const { return deforest_flag_; }
    
private:
    int deforest_flag_;  ///< Controls deForest prescription (0, 1, or -1)
};

/**
 * @class PionCrossSection
 * @brief Pion electroproduction cross section
 * 
 * Supports both exclusive (p(e,e'pi)X) and semi-inclusive production
 */
class PionCrossSection : public CrossSectionBase {
public:
    explicit PionCrossSection(PionType type = PionType::PI_PLUS);
    
    double Calculate(const SimcEvent& evt) const override;
    ReactionType GetReactionType() const override { return ReactionType::PION_PRODUCTION; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    /**
     * @brief Set pion type
     */
    void SetPionType(PionType type) { pion_type_ = type; }
    
private:
    PionType pion_type_;
    
    /**
     * @brief Get hadronic tensor for pion production
     */
    double HadronicTensor(const SimcEvent& evt) const;
};

/**
 * @class KaonCrossSection
 * @brief Kaon electroproduction cross section
 */
class KaonCrossSection : public CrossSectionBase {
public:
    explicit KaonCrossSection(KaonType type = KaonType::K_PLUS);
    
    double Calculate(const SimcEvent& evt) const override;
    ReactionType GetReactionType() const override { return ReactionType::KAON_PRODUCTION; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    void SetKaonType(KaonType type) { kaon_type_ = type; }
    
private:
    KaonType kaon_type_;
};

/**
 * @class CrossSectionFactory
 * @brief Factory to create appropriate cross section calculator
 */
class CrossSectionFactory {
public:
    /**
     * @brief Create cross section calculator based on reaction type
     * @param type Reaction type
     * @return Unique pointer to cross section calculator
     */
    static std::unique_ptr<CrossSectionBase> Create(ReactionType type);
    
    /**
     * @brief Create from configuration
     * @param config ConfigManager with reaction settings
     * @return Unique pointer to cross section calculator
     */
    static std::unique_ptr<CrossSectionBase> CreateFromConfig(const class ConfigManager& config);
};

/**
 * @class RadiativeCorrections
 * @brief Calculate radiative corrections to cross sections
 * 
 * Implements peaking approximation and other RC methods
 */
class RadiativeCorrections {
public:
    /**
     * @brief Constructor
     * @param bt Radiation thickness parameter (default 2/3)
     */
    explicit RadiativeCorrections(double bt = 0.667);
    
    /**
     * @brief Calculate radiative correction factor
     * @param Ein Incident energy (MeV)
     * @param Ee Scattered energy (MeV)
     * @param theta Scattering angle (rad)
     * @param Z Target atomic number
     * @return Correction factor (multiply cross section by this)
     */
    double GetCorrectionFactor(double Ein, double Ee, double theta, int Z) const;
    
    /**
     * @brief Calculate internal bremsstrahlung correction
     */
    double InternalCorrection(double Ein, double Ee, double theta) const;
    
    /**
     * @brief Calculate external bremsstrahlung (target)
     */
    double ExternalCorrection(double Ein, double Ee, int Z, double thickness) const;
    
private:
    double bt_;  // Radiation thickness parameter
    
    /**
     * @brief Calculate radiative tail
     */
    double RadiativeTail(double Ein, double Ee, double theta) const;
};

} // namespace simc

#endif // SIMC_CROSS_SECTION_H
