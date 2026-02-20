// tests/test_sos_matrix_debug.cpp
// Debug: Check if SOS matrix terms are actually loaded

#include "simc/spectrometers/SOS.h"
#include <iostream>
#include <iomanip>

using namespace simc;

int main() {
    std::cout << "=== SOS MATRIX LOADING DEBUG ===" << std::endl;
    std::cout << std::endl;
    
    SOS sos;
    
    if (!sos.LoadMatrices("../data/matrices/sos/forward_cosy.dat",
                          "../data/matrices/sos/recon_cosy.dat")) {
        std::cerr << "Failed to load matrices!" << std::endl;
        return 1;
    }
    
    std::cout << "Now testing forward transport with PURE Y angle:" << std::endl;
    std::cout << std::endl;
    
    SOS::TrackState track;
    track.x = 0.0;       // No X
    track.y = 0.0;       // Start at Y=0
    track.z = 0.0;
    track.dx = 0.0;      // No X angle
    track.dy = 0.010;    // 10 mrad Y angle (0.010 slope)
    track.delta = 0.0;   // No momentum deviation
    track.p = 1.0;
    track.m2 = 0.000511 * 0.000511;
    
    std::cout << "BEFORE transport:" << std::endl;
    std::cout << "  x=" << track.x << " cm, dx=" << track.dx << std::endl;
    std::cout << "  y=" << track.y << " cm, dy=" << track.dy << std::endl;
    std::cout << std::endl;
    
    // Transport through SOS
    bool success = sos.Transport(track);
    
    std::cout << "AFTER transport (success=" << (success ? "true" : "false") << "):" << std::endl;
    std::cout << "  Focal plane x=" << sos.GetFocalPlaneX() << " cm" << std::endl;
    std::cout << "  Focal plane dx=" << sos.GetFocalPlaneDX() << std::endl;
    std::cout << "  Focal plane y=" << sos.GetFocalPlaneY() << " cm" << std::endl;
    std::cout << "  Focal plane dy=" << sos.GetFocalPlaneDY() << std::endl;
    std::cout << std::endl;
    
    if (std::abs(sos.GetFocalPlaneY()) < 0.001) {
        std::cout << "ERROR: Y is still zero!" << std::endl;
        std::cout << std::endl;
        std::cout << "This means one of:" << std::endl;
        std::cout << "  1. Matrix terms are not loaded" << std::endl;
        std::cout << "  2. All transformations are marked as 'drift'" << std::endl;
        std::cout << "  3. Y coefficient terms are all zero" << std::endl;
        std::cout << "  4. TranspMatrix is not being called" << std::endl;
    } else {
        std::cout << "Good! Y changed to " << sos.GetFocalPlaneY() << " cm" << std::endl;
        std::cout << "For 10 mrad input, expect Y ~ 10-30 cm at focal plane" << std::endl;
    }
    
    return 0;
}
