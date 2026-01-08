// include/simc/CrossSection.h
// Cross section calculations for electron scattering

#ifndef SIMC_CROSS_SECTION_H
#define SIMC_CROSS_SECTION_H

#include "SimcTypes.h"
#include "SimcEvent.h"
#include <memory>

namespace simc {

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
 * Uses Rosenbluth formula with form factors
 */
class ElasticCrossSection : public CrossSectionBase {
public:
    ElasticCrossSection() = default;
    
    double Calculate(const SimcEvent& evt) const override;
    ReactionType GetReactionType() const override { return ReactionType::ELASTIC; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    /**
     * @brief Calculate Mott cross section
     * @param Ein Incident energy (MeV)
     * @param theta Scattering angle (rad)
     * @return d\sigma/d\Omega in ub/sr
     */
    static double MottCrossSection(double Ein, double theta);
    
    /**
     * @brief Calculate recoil factor
     * @param Ein Incident energy (MeV)
     * @param theta Scattering angle (rad)
     * @param M_target Target mass (MeV)
     * @return Recoil factor
     */
    static double RecoilFactor(double Ein, double theta, double M_target);
    
    /**
     * @brief Get electric form factor GE
     * @param Q2 Four-momentum transfer squared (GeV^2)
     * @return GE (dimensionless)
     */
    static double GetGE(double Q2);
    
    /**
     * @brief Get magnetic form factor GM
     * @param Q2 Four-momentum transfer squared (GeV^2)
     * @return GM (dimensionless)
     */
    static double GetGM(double Q2);
    
    /**
     * @brief Dipole form factor
     * @param Q2 in GeV^2
     * @param Lambda in GeV
     * @return Form factor value
     */
    static double DipoleFF(double Q2, double Lambda = 0.71);
};

/**
 * @class QuasiElasticCrossSection
 * @brief Quasi-elastic (e,e'p) scattering from nuclei
 * 
 * Includes Fermi motion, binding energy, and spectral function
 */
class QuasiElasticCrossSection : public CrossSectionBase {
public:
    QuasiElasticCrossSection() = default;
    
    double Calculate(const SimcEvent& evt) const override;
    ReactionType GetReactionType() const override { return ReactionType::QUASIELASTIC; }
    bool IsPhysical(const SimcEvent& evt) const override;
    
    /**
     * @brief Get off-shell correction factor
     * @param Pm Missing momentum (MeV/c)
     * @param Em Missing energy (MeV)
     * @return Correction factor
     */
    static double OffShellCorrection(double Pm, double Em);
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
