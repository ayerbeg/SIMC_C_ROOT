#include <iostream>
#include <iomanip>
#include <cmath>

// Test the angle conversions with EXACT Fortran formulas
// to find where the bug is

const double pi = 3.141592653589793;
const double deg_to_rad = pi / 180.0;

void physics_angles(double theta0, double phi0, double xptar, double yptar,
                    double& theta, double& phi) {
    // EXACT Fortran implementation from event.f lines 1669-1709
    
    double costh = std::cos(theta0);
    double sinth = std::sin(theta0);
    double cosph = std::cos(phi0);
    double sinph = std::sin(phi0);
    
    // Direction in spectrometer frame
    double dx = xptar;
    double dy = yptar;
    double dz = std::sqrt(1.0 + dx*dx + dy*dy);
    
    // Rotate to lab frame - EXACT Fortran formulas
    double ux = (costh * cosph * dx - sinph * dy + sinth * cosph * dz) / dz;
    double uy = (costh * sinph * dx + cosph * dy + sinth * sinph * dz) / dz;
    double uz = (-sinth * dx + costh * dz) / dz;
    
    std::cout << "  Direction in lab frame:" << std::endl;
    std::cout << "    ux = " << ux << std::endl;
    std::cout << "    uy = " << uy << std::endl;
    std::cout << "    uz = " << uz << std::endl;
    std::cout << "    |u| = " << std::sqrt(ux*ux + uy*uy + uz*uz) << std::endl;
    
    // Convert to angles
    theta = std::acos(uz);
    
    // Fortran lines 1703-1707:
    // if (dx.ne.0.0) then
    //   phi = atan( (dy*costh + sinth*sinph) / dx )
    //   if (phi.le.0) phi=phi+pi
    //   if (sinph.lt.0.) phi=phi+pi
    // else
    //   phi = phi0
    // endif
    
    if (dx != 0.0) {
        phi = std::atan((dy*costh + sinth*sinph) / dx);
        std::cout << "  Raw phi from atan = " << phi * 180.0/pi << " deg" << std::endl;
        
        if (phi <= 0.0) {
            phi = phi + pi;
            std::cout << "  After phi<=0 correction: " << phi * 180.0/pi << " deg" << std::endl;
        }
        
        if (sinph < 0.0) {
            phi = phi + pi;
            std::cout << "  After sinph<0 correction: " << phi * 180.0/pi << " deg" << std::endl;
        }
    } else {
        phi = phi0;
    }
}

void spectrometer_angles(double theta0, double phi0, double theta, double phi,
                        double& xptar, double& yptar) {
    // EXACT Fortran implementation from event.f lines 1720-1750
    
    double costh = std::cos(theta0);
    double sinth = std::sin(theta0);
    double cosph = std::cos(phi0);
    double sinph = std::sin(phi0);
    
    // Unit vector in lab frame
    double x = std::sin(theta) * std::cos(phi);
    double y = std::sin(theta) * std::sin(phi);
    double z = std::cos(theta);
    
    std::cout << "  Lab frame direction:" << std::endl;
    std::cout << "    x = " << x << std::endl;
    std::cout << "    y = " << y << std::endl;
    std::cout << "    z = " << z << std::endl;
    
    // Spectrometer central ray
    double x0 = sinth * cosph;
    double y0 = sinth * sinph;
    double z0 = costh;
    
    // cos(angle between event and central ray)
    double cos_dtheta = x*x0 + y*y0 + z*z0;
    
    std::cout << "  cos_dtheta = " << cos_dtheta << std::endl;
    std::cout << "  dtheta = " << std::acos(cos_dtheta) * 180.0/pi << " deg" << std::endl;
    
    // Fortran lines 1739-1740:
    // dx = x / cos_dtheta
    // dy = sqrt(1/cos_dtheta**2-1.-dx**2)
    xptar = x / cos_dtheta;
    
    double dy_squared = 1.0/(cos_dtheta*cos_dtheta) - 1.0 - xptar*xptar;
    std::cout << "  dy_squared = " << dy_squared << std::endl;
    
    if (dy_squared < 0.0) {
        std::cout << "  ERROR: dy_squared < 0!" << std::endl;
        dy_squared = 0.0;
    }
    
    yptar = std::sqrt(dy_squared);
    
    // Determine sign - Fortran lines 1745-1747:
    // y_event = y/cos_dtheta
    // if (y_event .lt. y0) dy = -dy
    double y_event = y / cos_dtheta;
    
    std::cout << "  y_event = " << y_event << std::endl;
    std::cout << "  y0 = " << y0 << std::endl;
    
    if (y_event < y0) {
        yptar = -yptar;
        std::cout << "  Negated yptar" << std::endl;
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    // Test case from your output
    double theta0 = 12.5 * deg_to_rad;
    double phi0 = 90.0 * deg_to_rad;
    
    double xptar_orig = 0.02;
    double yptar_orig = 0.01;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Angle Conversion Round-Trip Test" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Central spectrometer angles:" << std::endl;
    std::cout << "  theta0 = " << theta0 * 180.0/pi << " deg" << std::endl;
    std::cout << "  phi0 = " << phi0 * 180.0/pi << " deg" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Original spectrometer angles:" << std::endl;
    std::cout << "  xptar = " << xptar_orig << " rad" << std::endl;
    std::cout << "  yptar = " << yptar_orig << " rad" << std::endl;
    std::cout << std::endl;
    
    // Step 1: Convert to physics angles
    std::cout << "STEP 1: Spectrometer -> Physics" << std::endl;
    std::cout << "=================================" << std::endl;
    double theta, phi;
    physics_angles(theta0, phi0, xptar_orig, yptar_orig, theta, phi);
    
    std::cout << "Physics angles:" << std::endl;
    std::cout << "  theta = " << theta * 180.0/pi << " deg" << std::endl;
    std::cout << "  phi = " << phi * 180.0/pi << " deg" << std::endl;
    std::cout << std::endl;
    
    // Step 2: Convert back to spectrometer angles
    std::cout << "STEP 2: Physics -> Spectrometer" << std::endl;
    std::cout << "=================================" << std::endl;
    double xptar_back, yptar_back;
    spectrometer_angles(theta0, phi0, theta, phi, xptar_back, yptar_back);
    
    std::cout << "Recovered spectrometer angles:" << std::endl;
    std::cout << "  xptar = " << xptar_back << " rad" << std::endl;
    std::cout << "  yptar = " << yptar_back << " rad" << std::endl;
    std::cout << std::endl;
    
    // Check errors
    std::cout << "ERRORS:" << std::endl;
    std::cout << "========" << std::endl;
    double xptar_error = xptar_back - xptar_orig;
    double yptar_error = yptar_back - yptar_orig;
    
    std::cout << "  xptar error = " << xptar_error << " rad";
    std::cout << " (" << xptar_error * 1000.0 << " mrad)" << std::endl;
    std::cout << "  yptar error = " << yptar_error << " rad";
    std::cout << " (" << yptar_error * 1000.0 << " mrad)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "Relative errors:" << std::endl;
    std::cout << "  xptar: " << std::abs(xptar_error / xptar_orig) * 100.0 << "%" << std::endl;
    std::cout << "  yptar: " << std::abs(yptar_error / yptar_orig) * 100.0 << "%" << std::endl;
    
    return 0;
}
