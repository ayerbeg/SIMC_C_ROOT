// include/simc/spectrometers/HRSl.h
// Hall A Left HRS Spectrometer - Phase 5c.6

#ifndef SIMC_HRSR_H
#define SIMC_HRSR_H

#include <array>
#include <vector>
#include <string>

namespace simc {

/**
 * @brief Hall A Left High Resolution Spectrometer (HRSl)
 * 
 * Based on Jefferson Lab Hall A HRS-R spectrometer
 * Fortran reference: hrsr/mc_hrsr.f, hrsr/apertures_hrsr.inc
 * 
 * Configuration:
 * - 3 Quadrupoles (Q1: r=15.0 cm, Q2/Q3: r=30.22 cm)
 * - 1 Dipole with complex rotated apertures
 * - Collimator at z=110 cm (3.145 × 6.090 cm entrance)
 * - Momentum range: typically 0.3-4.0 GeV/c
 * - Momentum acceptance: ±4.5% (δp/p)
 * - 12 transformation classes
 */
class HRSl {
public:
    /**
     * @brief Track state during transport
     */
    struct TrackState {
        double x{0.0};        // Position X (cm)
        double y{0.0};        // Position Y (cm)
        double z{0.0};        // Position Z (cm)
        double dx{0.0};       // Slope dX/dZ
        double dy{0.0};       // Slope dY/dZ
        double delta{0.0};    // Momentum deviation (%)
        double p{0.0};        // Momentum (GeV/c)
        double m2{0.0};       // Mass squared (GeV²)
        double pathlen{0.0};  // Path length traveled (cm)
        bool decayed{false};  // Has particle decayed?
    };

    struct FocalPlaneState {
        double x{0.0};      // X at focal plane (cm)
        double y{0.0};      // Y at focal plane (cm)  
        double dx{0.0};     // dX/dZ at focal plane
        double dy{0.0};     // dY/dZ at focal plane
    };

    struct TargetState {
        double xptar{0.0};  // X angle at target (rad)
        double yptar{0.0};  // Y angle at target (rad)
        double delta{0.0};  // Momentum deviation (%)
        double ytar{0.0};   // Y position at target (cm)
    };

    /**
     * @brief Aperture structures
     */
    struct CircularAperture {
        double radius{0.0};  // Aperture radius (cm)
    };

    struct RectangularAperture {
        double x_min{0.0};
        double x_max{0.0};
        double y_min{0.0};
        double y_max{0.0};
    };

    /**
     * @brief Matrix element structure
     */
    struct MatrixTerm {
        double coeff[5]{0.0};
        int expon[5]{0};
    };

    struct Transformation {
        std::vector<MatrixTerm> terms;
        int n_terms{0};
        double length{0.0};
        bool is_drift{false};
    };

    struct MatrixElements {
        static constexpr int MAX_CLASSES = 20;
        static constexpr int MAX_COEFFS = 1000;
        std::array<Transformation, MAX_CLASSES> classes;
        int n_classes{0};
    };

    /**
     * @brief Transport statistics (matches Fortran lSTOP_* counters)
     */
    struct TransportStats {
        int lSTOP_trials{0};
        int lSTOP_successes{0};
        int lSTOP_slit{0};         // Collimator rejections
        int lSTOP_Q1_in{0};
        int lSTOP_Q1_mid{0};
        int lSTOP_Q1_out{0};
        int lSTOP_Q2_in{0};
        int lSTOP_Q2_mid{0};
        int lSTOP_Q2_out{0};
        int lSTOP_D1_in{0};
        int lSTOP_D1_out{0};
        int lSTOP_Q3_in{0};
        int lSTOP_Q3_mid{0};
        int lSTOP_Q3_out{0};
        int lSTOP_hut{0};

        void Reset() {
            lSTOP_trials = 0;
            lSTOP_successes = 0;
            lSTOP_slit = 0;
            lSTOP_Q1_in = 0;
            lSTOP_Q1_mid = 0;
            lSTOP_Q1_out = 0;
            lSTOP_Q2_in = 0;
            lSTOP_Q2_mid = 0;
            lSTOP_Q2_out = 0;
            lSTOP_D1_in = 0;
            lSTOP_D1_out = 0;
            lSTOP_Q3_in = 0;
            lSTOP_Q3_mid = 0;
            lSTOP_Q3_out = 0;
            lSTOP_hut = 0;
        }
    };

    // Constructor
    HRSl();

    // Main transport and reconstruction
    bool Transport(TrackState& track);
    bool Reconstruct(const FocalPlaneState& fp, TargetState& tgt);

    // Matrix file loading
    bool LoadMatrices(const std::string& forward_file, const std::string& recon_file);

    // Statistics
    const TransportStats& GetStats() const { return stats_; }
    void ResetStats() { stats_.Reset(); }

private:
    // Apertures (from apertures_hrsr.inc and mc_hrsr.f)
    CircularAperture q1_;           // Q1: r=15.0 cm
    CircularAperture q2_;           // Q2: r=30.22 cm  
    CircularAperture q3_;           // Q3: r=30.22 cm

    // Statistics
    TransportStats stats_;

    // Matrix elements
    MatrixElements forward_matrix_;
    MatrixElements recon_matrix_;

    // Transport helper functions  
    void TranspMatrix(TrackState& track, int class_num, double pathlen);
    void Project(TrackState& track, double distance);
    
    // Aperture checks
    bool CheckCircularAperture(double x, double y, double r_sq);
    bool CheckRectangularAperture(double x, double y, const RectangularAperture& ap);
    
    // Matrix parsing
    bool ParseMatrixFile(const std::string& filename, MatrixElements& matrix);
};

} // namespace simc

#endif // SIMC_HRSR_H
