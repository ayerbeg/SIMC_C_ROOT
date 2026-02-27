#ifndef SIMC_SPECTROMETERS_SHMS_H
#define SIMC_SPECTROMETERS_SHMS_H

#include <string>
#include <vector>
#include <array>

namespace simc {

// Forward declarations
struct TrackState;
struct FocalPlaneState;
struct TargetState;

/**
 * @brief SHMS (Super High Momentum Spectrometer) Transport and Reconstruction
 * 
 * Implements Monte Carlo transport through the SHMS spectrometer including:
 * - 31 aperture checks (4 HB, 15 Quad, 12 Dipole)
 * - Matrix element transport through magnetic elements
 * - Field-free drift calculations
 * - Coordinate rotations for tilted apertures
 * - Focal plane reconstruction
 * 
 * Based on simc_gfortran mc_shms.f (2017 ME version)
 * 
 * **PHASE 5e MODIFICATIONS:**
 * - Added GetFocalPlane() method to extract focal plane coordinates after transport
 */
class SHMS {
public:
    /**
     * @brief Track state during transport
     */
    struct TrackState {
        double x{0.0};        // Position X (cm)
        double y{0.0};        // Position Y (cm)
        double z{0.0};        // Position Z (cm)
        double dx{0.0};       // Slope dX/dZ (dimensionless)
        double dy{0.0};       // Slope dY/dZ (dimensionless)
        double delta{0.0};    // Momentum deviation (%)
        double p{0.0};        // Momentum (GeV/c)
        double m2{0.0};       // Mass squared (GeV²)
        double pathlen{0.0};  // Path length traveled (cm)
        bool decayed{false};  // Has particle decayed?
    };

    /**
     * @brief Focal plane coordinates
     */
    struct FocalPlaneState {
        double x{0.0};      // X at focal plane (cm)
        double xp{0.0};     // X' = dX/dZ (mrad)
        double y{0.0};      // Y at focal plane (cm)
        double yp{0.0};     // Y' = dY/dZ (mrad)
        double delta{0.0};  // Momentum deviation (%)
    };

    /**
     * @brief Target coordinates (reconstructed)
     */
    struct TargetState {
        double x{0.0};      // X at target (cm)
        double xp{0.0};     // X' = dX/dZ (mrad)
        double y{0.0};      // Y at target (cm)
        double yp{0.0};     // Y' = dY/dZ (mrad)
        double delta{0.0};  // Momentum deviation (%)
    };

    /**
     * @brief Drift lengths for SHMS (2017 ME version, 26cm offset)
     */
    struct DriftLengths {
        // Holding Box (Bender)
        double hb_in{118.39};      // Target to HB entrance
        double hb_men{17.61};      // HB entrance to mag entrance
        double hb_mex{80.00};      // HB mag entrance to mag exit (MATRIX)
        double hb_out{17.61};      // HB mag exit to exit
        
        // Q1 (Quadrupole 1)
        double q1_in{58.39};       // HB exit to Q1 entrance
        double q1_men{28.35};      // Q1 entrance to mag entrance
        double q1_mid{93.65};      // Q1 mag entrance to mid (MATRIX)
        double q1_mex{93.65};      // Q1 mid to mag exit (MATRIX)
        double q1_out{28.35};      // Q1 mag exit to exit
        
        // Q2 (Quadrupole 2)
        double q2_in{25.55};       // Q1 exit to Q2 entrance
        double q2_men{39.10};      // Q2 entrance to mag entrance
        double q2_mid{79.35};      // Q2 mag entrance to mid (MATRIX)
        double q2_mex{79.35};      // Q2 mid to mag exit (MATRIX)
        double q2_out{39.10};      // Q2 mag exit to exit
        
        // Q3 (Quadrupole 3)
        double q3_in{28.10};       // Q2 exit to Q3 entrance
        double q3_men{39.10};      // Q3 entrance to mag entrance
        double q3_mid{79.35};      // Q3 mag entrance to mid (MATRIX)
        double q3_mex{79.35};      // Q3 mid to mag exit (MATRIX)
        double q3_out{39.10};      // Q3 mag exit to exit
        
        // Dipole (D1)
        double q3_d1_trans{18.00}; // Q3 exit to D1 transition
        double d1_flare{30.10};    // D1 flare section
        double d1_men{39.47};      // D1 mag entrance (MATRIX)
        double d1_mid[7];          // D1 mid sections (MATRIX, each 36.406263)
        double d1_mex{36.406263};  // D1 mag exit (MATRIX)
        double d1_out{60.68};      // D1 exit
        
        // Focal Plane
        double fp{307.95};         // Final drift to focal plane
        
        DriftLengths() {
            // Initialize D1 mid sections
            for (int i = 0; i < 7; ++i) {
                d1_mid[i] = 36.406263;
            }
        }
    };

    /**
     * @brief Aperture definitions for SHMS
     */
    struct Apertures {
        // Holding Box (asymmetric, tilted)
        double r_HBx{11.2};                              // Vertical half-width (cm)
        std::array<double, 4> r_HBfyp{11.75, 11.74, 11.71, 11.70};  // Upper y limits
        std::array<double, 4> r_HBfym{-4.13, -5.45, -10.25, -11.71}; // Lower y limits
        std::array<double, 4> hb_tilt{1.5, 1.5, -1.5, -1.5};        // Tilt angles (deg)
        std::array<double, 4> hb_offset{1.51, 0.98, 0.98, 1.51};    // Y offsets (cm)
        
        // Quadrupoles (circular)
        double r_Q1{20.0};  // Q1 radius (cm)
        double r_Q2{30.0};  // Q2 radius (cm)
        double r_Q3{30.0};  // Q3 radius (cm)
        
        // Dipole (circular, tilted, offset)
        double r_D1{30.0};  // D1 radius (cm)
        std::array<double, 12> d1_angles{
            9.2, 6.9, 4.6, 2.3, 0.0, -2.3, -4.6, -6.9, -9.2, 0.0, 0.0, 0.0
        };  // Tilt angles (deg)
        std::array<double, 12> d1_offsets{
            14.70, 11.02, 7.35, 3.67, 0.0, -3.67, -7.35, -11.02, -14.70, 0.0, 0.0, 0.0
        };  // X offsets (cm)
        
        // Collimator (optional, octagonal)
        bool use_collimator{false};
        double coll_x_in{8.5};    // Entrance half-width (cm)
        double coll_y_in{12.5};   // Entrance half-height (cm)
        double coll_x_out{8.65};  // Exit half-width (cm)
        double coll_y_out{12.85}; // Exit half-height (cm)
    };

    /**
     * @brief Matrix element term
     */
    struct MatrixTerm {
        std::array<double, 5> coeff{0.0};  // Coefficients for X, XP, Y, YP, dL
        std::array<int, 5> expon{0};       // Exponents for x, xp, y, yp, delta
    };

    /**
     * @brief Matrix transformation class
     */
    struct Transformation {
        int n_terms{0};
        std::vector<MatrixTerm> terms;
        double length{0.0};  // Canonical length (cm)
    };

    /**
     * @brief Matrix elements for transport
     */
    struct MatrixElements {
        static constexpr int MAX_CLASSES = 32;  // SHMS has 32 transformation classes
        std::array<Transformation, MAX_CLASSES> classes;
    };

    /**
     * @brief Transport statistics (aperture rejection counters)
     * Member names match Fortran SIMC convention (shmsSTOP_*)
     */
    struct TransportStats {
        int total_events{0};
        int accepted{0};
        
        // Aggregate counters (sum of all apertures in each element)
        int shmsSTOP_HB{0};     // Total rejected by Holding Box
        int shmsSTOP_Q1{0};     // Total rejected by Q1
        int shmsSTOP_Q2{0};     // Total rejected by Q2
        int shmsSTOP_Q3{0};     // Total rejected by Q3
        int shmsSTOP_D1{0};     // Total rejected by Dipole
        int shmsSTOP_COLL{0};   // Total rejected by Collimator
        int shmsSTOP_HUT{0};    // Total rejected in Hut
        
        // Detailed counters (individual apertures)
        // Holding Box
        int stop_hb_in{0};
        int stop_hb_men{0};
        int stop_hb_mex{0};
        int stop_hb_out{0};
        
        // Q1
        int stop_q1_in{0};
        int stop_q1_men{0};
        int stop_q1_mid{0};
        int stop_q1_mex{0};
        int stop_q1_out{0};
        
        // Q2
        int stop_q2_in{0};
        int stop_q2_men{0};
        int stop_q2_mid{0};
        int stop_q2_mex{0};
        int stop_q2_out{0};
        
        // Q3
        int stop_q3_in{0};
        int stop_q3_men{0};
        int stop_q3_mid{0};
        int stop_q3_mex{0};
        int stop_q3_out{0};
        
        // Dipole
        int stop_d1_in{0};     // D1 mechanical entrance (Q3/D1 transition)
        int stop_d1_flare{0};
        int stop_d1_men{0};
        int stop_d1_mid1{0};
        int stop_d1_mid2{0};
        int stop_d1_mid3{0};
        int stop_d1_mid4{0};
        int stop_d1_mid5{0};
        int stop_d1_mid6{0};
        int stop_d1_mid7{0};
        int stop_d1_mex{0};
        int stop_d1_out{0};
        
        // Collimator
        int stop_coll_in{0};
        int stop_coll_out{0};
        
        void Reset() {
            *this = TransportStats();
        }
        
        double GetAcceptance() const {
            return total_events > 0 ? 
                   static_cast<double>(accepted) / total_events : 0.0;
        }
    };

    // Constructor/Destructor
    SHMS();
    ~SHMS() = default;

    // Main interface
    /**
     * @brief Transport particle through SHMS spectrometer
     * @param track Track state (updated in place)
     * @return true if particle accepted, false if rejected by aperture
     */
    bool Transport(TrackState& track);

    /**
     * @brief Reconstruct target coordinates from focal plane
     * @param fp Focal plane state
     * @param target Target state (output)
     * @return true if reconstruction successful
     */
    bool Reconstruct(const FocalPlaneState& fp, TargetState& target);

    // ========================================================================
    // PHASE 5e: FOCAL PLANE EXTRACTION
    // ========================================================================
    /**
     * @brief Extract focal plane coordinates after transport
     * @param track Track state at focal plane (after successful Transport() call)
     * @param fp Focal plane state to fill (output)
     * @return true if extraction successful
     * 
     * This method should be called AFTER Transport() succeeds.
     * It extracts the focal plane coordinates from the transported track state.
     * 
     * UNITS:
     * - x, y: cm
     * - xp, yp: mrad (converted from slopes)
     * - delta: %
     * 
     * Phase 5e - Week 1, Day 1
     */
    bool GetFocalPlane(const TrackState& track, FocalPlaneState& fp) const;

    // Configuration
    /**
     * @brief Load matrix files for forward transport and reconstruction
     * @param forward_file Path to forward matrix file
     * @param recon_file Path to reconstruction matrix file
     * @return true if successful
     */
    bool LoadMatrices(const std::string& forward_file, const std::string& recon_file);

    /**
     * @brief Set whether to use collimator
     */
    void SetUseCollimator(bool use) { apertures_.use_collimator = use; }

    // Statistics
    const TransportStats& GetStatistics() const { return stats_; }
    const TransportStats& GetStats() const { return stats_; }  // Alias for convenience
    void ResetStatistics() { stats_.Reset(); }

private:
    // Element transport functions
    bool TransportHoldingBox(TrackState& track);
    bool TransportQ1(TrackState& track);
    bool TransportQ2(TrackState& track);
    bool TransportQ3(TrackState& track);
    bool TransportDipole(TrackState& track);
    bool TransportHut(TrackState& track);

    // Aperture checking
    bool CheckHBAperture(double x, double y, int aperture_id);
    bool CheckQuadAperture(double x, double y, double radius);
    bool CheckDipoleAperture(double x, double y, int section);
    bool CheckCollimatorAperture(double x, double y, bool entrance);

    // Helper functions
    void Project(TrackState& track, double zdrift);
    void TranspMatrix(TrackState& track, int class_id, double zdrift);
    void RotateVAxis(double angle_deg, double& x, double& y, double dx, double dy);
    void RotateHAxis(double angle_deg, double& x, double& y, double dx, double dy);

    // Matrix file parsing
    bool ParseMatrixFile(const std::string& filepath, MatrixElements& matrices);

    // Data members
    DriftLengths drifts_;
    Apertures apertures_;
    MatrixElements forward_matrix_;
    MatrixElements recon_matrix_;
    TransportStats stats_;
};

} // namespace simc

#endif // SIMC_SPECTROMETERS_SHMS_H
