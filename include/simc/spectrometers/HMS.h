#ifndef SIMC_SPECTROMETERS_HMS_H
#define SIMC_SPECTROMETERS_HMS_H

#include <string>
#include <vector>
#include <array>

namespace simc {

// Forward declarations
struct TrackState;
struct FocalPlaneState;
struct TargetState;

/**
 * @brief HMS (High Momentum Spectrometer) Transport and Reconstruction
 * 
 * Implements Monte Carlo transport through the HMS spectrometer including:
 * - Collimator (octagonal, optional)
 * - 3 Quadrupoles (Q1, Q2, Q3)
 * - 1 Dipole with complex shape
 * - Vacuum pipes post-dipole
 * - Matrix element transport through magnetic elements
 * - Focal plane reconstruction
 * 
 * Based on simc_gfortran hms/mc_hms.f (HMS-100 tune)
 * 
 * PHASE 5e: Added GetFocalPlane() method - Week 1, Day 1
 */
class HMS {
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
     * @brief HMS-100 Collimator dimensions (octagonal)
     * Can be disabled for testing
     */
    struct CollimatorDimensions {
        bool use_collimator{true};     // Enable/disable collimator
        double h_entr{4.575};          // Horizontal half-width at entrance (cm)
        double v_entr{11.646};         // Vertical half-width at entrance (cm)
        double h_exit{4.759};          // Horizontal half-width at exit (cm)
        double v_exit{12.114};         // Vertical half-width at exit (cm)
        double x_off{0.000};           // X offset (cm), +ve = DOWN
        double y_off{0.028};           // Y offset (cm), +ve = LEFT
        double z_off{40.17};           // Z offset for HMS-100 tune (cm)
        double z_entr{126.2 + 40.17};  // Z position of entrance (cm)
        double z_exit{126.2 + 40.17 + 6.3};  // Z position of exit (6.3 cm thick)
    };

    /**
     * @brief HMS apertures (quadrupoles and dipole)
     */
    struct Apertures {
        // Quadrupoles (circular)
        double r_Q1{20.50};  // Q1 radius (cm)
        double r_Q2{30.22};  // Q2 radius (cm)
        double r_Q3{30.22};  // Q3 radius (cm)
        
        // Dipole (complex shape - Niculescu model)
        // Dipole aperture defined as piecewise regions
        // See hit_dipole() function for full geometry
        double x_d1{34.29};   // Region 1 X limit (cm)
        double y_d1{12.07};   // Region 1 Y limit (cm)
        double x_d2{27.94};   // Region 2 X limit (cm)
        double y_d2{18.42};   // Region 2 Y limit (cm)
        double x_d3{13.97};   // Region 3 X limit (cm)
        double y_d3{18.95};   // Region 3 Y limit (cm)
        double x_d4{1.956};   // Region 4 X limit (cm)
        double y_d4{20.32};   // Region 4 Y limit (cm)
        double x_d5{27.94};   // Rounded corner center X (cm)
        double y_d5{12.065};  // Rounded corner center Y (cm)
        double r_d5{6.35};    // Rounded corner radius (cm)
        double a_d6{-0.114};  // Slant line slope
        double b_d6{20.54};   // Slant line Y-intercept (cm)
        
        // Post-dipole vacuum pipes
        double pipe_x_offset{2.8};   // X offset for pipes (cm)
        double pipe_y_offset{0.0};   // Y offset for pipes (cm)
        
        // Pipe positions (from dipole exit)
        double z_pipe1{64.77};       // End of 26.65" pipe (cm)
        double z_pipe2{64.77 + 297.18};  // End of 117" pipe (cm)
        double z_pipe3{64.77 + 297.18 + 115.57};  // End of 45.5" pipe (cm)
        
        // Pipe radii squared (for fast checking)
        double r_pipe1_sq{1145.518};    // (33.84 cm)² at exit of 26.65" pipe
        double r_pipe2_sq{1512.2299};   // (38.89 cm)² at exit of 117" pipe
        double r_pipe3_sq{2162.9383};   // (46.51 cm)² at exit of 45.5" pipe
    };

    /**
     * @brief HMS drift lengths
     * Based on HMS-100 tune with matrix transformations
     */
    struct DriftLengths {
        // Transformation 1: Target to Q1 entrance (drift)
        double target_to_q1;  // Read from matrix file
        
        // Transformation 2-3: Q1 (through quadrupole via matrix)
        double q1_mid{125.233};  // Q1 2/3 point length (cm)
        double q1_exit{62.617};  // Q1 mid to exit length (cm)
        
        // Transformation 4: Q1 to Q2 (drift)
        double q1_to_q2;     // Read from matrix file
        
        // Transformation 5-6: Q2 (through quadrupole via matrix)
        double q2_mid{143.90};  // Q2 2/3 point length (cm)
        double q2_exit{71.95};  // Q2 mid to exit length (cm)
        
        // Transformation 7: Q2 to Q3 (drift)
        double q2_to_q3;     // Read from matrix file
        
        // Transformation 8-9: Q3 (through quadrupole via matrix)
        double q3_mid{143.8};   // Q3 2/3 point length (cm)
        double q3_exit{71.9};   // Q3 mid to exit length (cm)
        
        // Transformation 10: Q3 to D1 (drift)
        double q3_to_d1;     // Read from matrix file
        
        // Transformation 11: D1 (through dipole via matrix)
        double d1_length{526.053};  // Dipole magnetic length (cm)
        
        // Transformation 12: D1 exit to focal plane (drift)
        double d1_to_fp;     // Read from matrix file
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
        double length{0.0};   // Canonical length (cm)
        bool is_drift{false}; // True if this is a drift space
    };

    /**
     * @brief Matrix elements for transport
     */
    struct MatrixElements {
        static constexpr int MAX_CLASSES = 20;  // HMS matrix files (increased for safety)
        std::array<Transformation, MAX_CLASSES> classes;
    };

    /**
     * @brief Transport statistics (aperture rejection counters)
     * Member names match Fortran SIMC convention (hSTOP_*)
     */
    struct TransportStats {
        int total_events{0};
        int accepted{0};
        
        // Collimator
        int hSTOP_slit_hor{0};   // Horizontal slit
        int hSTOP_slit_vert{0};  // Vertical slit
        int hSTOP_slit_oct{0};   // Octagonal corners
        int hSTOP_coll{0};       // Collimator total (if using HMS collimator simulation)
        
        // Quadrupoles
        int hSTOP_Q1_in{0};
        int hSTOP_Q1_mid{0};
        int hSTOP_Q1_out{0};
        int hSTOP_Q2_in{0};
        int hSTOP_Q2_mid{0};
        int hSTOP_Q2_out{0};
        int hSTOP_Q3_in{0};
        int hSTOP_Q3_mid{0};
        int hSTOP_Q3_out{0};
        
        // Dipole
        int hSTOP_D1_in{0};
        int hSTOP_D1_out{0};
        
        // Detector hut
        int hSTOP_hut{0};
        int hSTOP_dc1{0};
        int hSTOP_dc2{0};
        int hSTOP_scin{0};
        int hSTOP_cal{0};
        
        // Successes
        int hSTOP_trials{0};
        int hSTOP_successes{0};
        
        void Reset() {
            *this = TransportStats();
        }
        
        double GetAcceptance() const {
            return hSTOP_trials > 0 ? 
                   static_cast<double>(hSTOP_successes) / hSTOP_trials : 0.0;
        }
    };

    // Constructor/Destructor
    HMS();
    ~HMS() = default;

    // Main interface
    /**
     * @brief Transport particle through HMS spectrometer
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
    // PHASE 5e: FOCAL PLANE EXTRACTION - Week 1, Day 1
    // ========================================================================
    /**
     * @brief Extract focal plane coordinates after transport
     * @param track Track state at focal plane (after successful Transport())
     * @param fp Focal plane state to fill (output)
     * @return true if extraction successful
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
    void SetUseCollimator(bool use) { collimator_.use_collimator = use; }

    /**
     * @brief Enable/disable HMS collimator simulation
     * @param use If true, use detailed pion absorption model (not implemented yet)
     */
    void SetUseHMSCollimator(bool use) { use_hms_coll_sim_ = use; }

    // Statistics
    const TransportStats& GetStatistics() const { return stats_; }
    const TransportStats& GetStats() const { return stats_; }  // Alias
    void ResetStatistics() { stats_.Reset(); }

private:
    // Element transport functions
    bool TransportCollimator(TrackState& track);
    bool TransportQ1(TrackState& track);
    bool TransportQ2(TrackState& track);
    bool TransportQ3(TrackState& track);
    bool TransportDipole(TrackState& track);
    bool TransportPipes(TrackState& track);
    bool TransportHut(TrackState& track);

    // Aperture checking
    bool CheckCollimatorEntrance(double x, double y);
    bool CheckCollimatorExit(double x, double y);
    bool CheckQuadAperture(double x, double y, double radius);
    bool CheckDipoleAperture(double x, double y);
    bool CheckPipeAperture(double x, double y, double x_offset, double y_offset, double r_sq);

    // Helper functions
    void Project(TrackState& track, double zdrift);
    void TranspMatrix(TrackState& track, int class_id, double zdrift);
    void RotateHAxis(double angle_deg, double& x, double& y);
    
    /**
     * @brief Check if particle hits dipole aperture (complex shape)
     * Based on Niculescu model with 6 regions
     * @param x X position (cm)
     * @param y Y position (cm)
     * @return true if OUTSIDE aperture (hit), false if inside (pass)
     */
    bool HitDipole(double x, double y);

    // Matrix file parsing
    bool ParseMatrixFile(const std::string& filepath, MatrixElements& matrices);

    // Data members
    CollimatorDimensions collimator_;
    Apertures apertures_;
    DriftLengths drifts_;
    MatrixElements forward_matrix_;
    MatrixElements recon_matrix_;
    TransportStats stats_;
    bool use_hms_coll_sim_{false};  // Use detailed HMS collimator simulation
};

} // namespace simc

#endif // SIMC_SPECTROMETERS_HMS_H
