// include/simc/spectrometers/SOS.h
// SOS (Short Orbit Spectrometer) - Phase 5c.5
// Simplest Hall C spectrometer: 1 Quad + 2 Dipoles

#ifndef SIMC_SPECTROMETERS_SOS_H
#define SIMC_SPECTROMETERS_SOS_H

#include <string>
#include <vector>
#include <array>

namespace simc {

/**
 * @brief SOS (Short Orbit Spectrometer) Transport and Reconstruction
 * 
 * Simplest Hall C spectrometer with large acceptance:
 * - 1 Quadrupole (entrance, mid, exit)
 * - 2 Bending Magnets (dipoles)
 * - Large momentum acceptance: ±40%
 * - Momentum range: 0.17-1.7 GeV/c
 * - Solid angle: ~7 msr
 * 
 * Based on simc_gfortran sos/mc_sos.f
 */
class SOS {
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
        double y{0.0};      // Y at focal plane (cm)
        double dx{0.0};     // dX/dZ at focal plane
        double dy{0.0};     // dY/dZ at focal plane
    };

    /**
     * @brief Target coordinates
     */
    struct TargetState {
        double xptar{0.0};  // X angle at target (rad)
        double yptar{0.0};  // Y angle at target (rad)
        double delta{0.0};  // Momentum deviation (%)
        double ytar{0.0};   // Y position at target (cm)
    };

    /**
     * @brief Quadrupole aperture dimensions
     */
    struct QuadAperture {
        double radius{0.0};  // Aperture radius (cm)
    };

    /**
     * @brief Dipole aperture dimensions
     */
    struct DipoleAperture {
        double x_min{0.0};
        double x_max{0.0};
        double y_min{0.0};
        double y_max{0.0};
    };

    /**
     * @brief Drift lengths between elements
     */
    struct DriftLengths {
        double target_to_quad{0.0};    // Drift from target to Q entrance
        double quad_entrance{0.0};      // Q entrance region
        double quad_mid{0.0};           // Q mid region  
        double quad_exit{0.0};          // Q exit region
        double quad_to_bm1{0.0};        // Drift from Q to BM1
        double bm1_length{0.0};         // BM1 transport length
        double bm1_to_bm2{0.0};         // Drift between dipoles
        double bm2_length{0.0};         // BM2 transport length
        double bm2_to_fp{0.0};          // Drift to focal plane
    };

    /**
     * @brief Matrix element structure (COSY format)
     */
    struct MatrixTerm {
        double coeff[5]{0.0};  // 5 coefficients
        int expon[5]{0};       // 5 exponents for x, dx, y, dy, delta
    };

    struct Transformation {
        std::vector<MatrixTerm> terms;
        int n_terms{0};
        double length{0.0};
        bool is_drift{false};
    };

    struct MatrixElements {
        static constexpr int MAX_CLASSES = 20;  // SOS has fewer than HMS
        static constexpr int MAX_COEFFS = 1000;
        std::array<Transformation, MAX_CLASSES> classes;
        int n_classes{0};
    };

    /**
     * @brief Transport statistics (matches Fortran sSTOP_* counters)
     */
    struct TransportStats {
        int sSTOP_trials{0};       // Total events attempted
        int sSTOP_successes{0};    // Events that made it through
        int sSTOP_slit{0};         // Stopped at collimator slit
        int sSTOP_quad_in{0};      // Stopped at quad entrance
        int sSTOP_quad_mid{0};     // Stopped at quad middle
        int sSTOP_quad_out{0};     // Stopped at quad exit
        int sSTOP_bm1_in{0};       // Stopped at BM1 entrance
        int sSTOP_bm1_out{0};      // Stopped at BM1 exit
        int sSTOP_bm2_in{0};       // Stopped at BM2 entrance
        int sSTOP_bm2_out{0};      // Stopped at BM2 exit
        int sSTOP_hut{0};          // Made it to detector hut
        
        void Reset() {
            sSTOP_trials = 0;
            sSTOP_successes = 0;
            sSTOP_slit = 0;
            sSTOP_quad_in = 0;
            sSTOP_quad_mid = 0;
            sSTOP_quad_out = 0;
            sSTOP_bm1_in = 0;
            sSTOP_bm1_out = 0;
            sSTOP_bm2_in = 0;
            sSTOP_bm2_out = 0;
            sSTOP_hut = 0;
        }
    };

    // Constructor
    SOS();

    // Main transport function
    bool Transport(TrackState& track);

    // Reconstruction
    void Reconstruct(const FocalPlaneState& fp, TargetState& target);


    /**
     * @brief Extract focal plane coordinates after transport
     * Phase 5e - Week 1, Day 1
     */
    bool GetFocalPlane(const TrackState& track, FocalPlaneState& fp) const;
  
    // Matrix loading
    bool LoadMatrices(const std::string& forward_file, 
                      const std::string& recon_file);

    // Statistics
    const TransportStats& GetStats() const { return stats_; }
    void ResetStats() { stats_.Reset(); }

private:
    // Element transport functions
    bool TransportQuadrupole(TrackState& track);
    bool TransportBM1(TrackState& track);
    bool TransportBM2(TrackState& track);
    bool TransportHut(TrackState& track);

    // Aperture checking
    bool CheckQuadAperture(double x, double y, double r_sq);
    bool CheckDipoleAperture(double x, double y, const DipoleAperture& aperture);

    // Helper functions
    void Project(TrackState& track, double drift_length);
    void TranspMatrix(TrackState& track, int class_index, double drift_length);

    // Matrix file parsing
    bool ParseMatrixFile(const std::string& filepath, MatrixElements& matrices);

    // Configuration
    QuadAperture quad_;
    DipoleAperture bm1_;
    DipoleAperture bm2_;
    DriftLengths drifts_;

    // Matrix elements
    MatrixElements forward_matrix_;
    MatrixElements recon_matrix_;

    // Statistics
    TransportStats stats_;
};

} // namespace simc

#endif // SIMC_SPECTROMETERS_SOS_H
