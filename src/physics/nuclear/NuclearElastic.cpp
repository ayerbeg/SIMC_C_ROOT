// src/physics/nuclear/NuclearElastic.cpp — see include/simc/physics/nuclear/NuclearElastic.h.
// A(e,e')A nuclear elastic scattering (task U.1): the cross section of upstream SIMC nuclear_elastic.f and the
// IPhysicsModel "nuc_elastic" used for doing_nuc_elast (event.f:complete_main l.1457-1459 at commit 20d46fc).
// Arithmetic as gfortran -fdefault-real-8: every constant is real*8; x**2 is x*x; A**(1./3.) is pow(A, 1/3).
#include "simc/physics/nuclear/NuclearElastic.h"

#include <cmath>
#include <string>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/physics/proton/ProtonCrossSections.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/StageContext.h"

namespace simc {
namespace nuclear {

// Fortran: nuclear_elastic.f:Fgauss (l.48-66)
double Fgauss(double T, double A) {
    using constants::hbarc;
    double fg = 0.;
    double radius = 1.07 * std::pow(A, 1. / 3.);
    // from H. de Vries, Atomic Data and Nuclear Data Tables 36, 495 (1987)
    const long iA = std::lround(A);
    if (iA == 205) radius = 5.470;
    if (iA == 56) radius = 3.729;
    if (iA == 28) radius = 3.085;
    if (iA == 27) radius = 3.035;
    const double hbarc2 = hbarc * hbarc;   // hbarc**2 (folded constant)
    const double x2 = (T / hbarc2) * (radius * radius);
    const double chr = (T / hbarc2) * (2.4 * 2.4) / 6.;
    if (chr < 80) fg = std::exp(-chr) / (1. + x2 / 6.);
    return fg;
}

// Fortran: nuclear_elastic.f:NUC_FORM_FACTOR (l.30-45)
void NucFormFactor(double QSQ, double& W1, double& W2, double& FF, double Z, double A) {
    W1 = 0.;
    FF = Fgauss(QSQ, A);
    W2 = (Z * FF) * (Z * FF);
}

// Fortran: nuclear_elastic.f:sig_nuc_elastic (l.1-27)
double SigNucElastic(double Ein, double Eprime, double theta, double Q2, double Z, double A, double M) {
    double W1, W2, FF;
    NucFormFactor(Q2, W1, W2, FF, Z, A);
    const double csmott = proton::sigMott(Eprime, theta, Q2);
    const double recoil = M / (M + Ein * (1. - std::cos(theta)));
    const double t = std::tan(theta / 2.);
    return (W2 + 2. * W1 * (t * t)) * csmott * recoil;
}

double SigNucElastic(const Event& ev, double Z, double A, double M) {
    return SigNucElastic(ev.Ein, ev.e.E, ev.e.theta, ev.Q2, Z, A, M);
}

namespace {
/// IPhysicsModel "nuc_elastic": event.f:complete_main l.1457-1459 (20d46fc), sigcc = sig_nuc_elastic(vertex/recon).
class NucElasticModel : public IPhysicsModel {
public:
    std::string name() const override { return "nuc_elastic"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext& ctx) override {
        Z_ = ctx.setup.targ.Z;
        A_ = ctx.setup.targ.A;
        M_ = ctx.setup.targ.M;
    }
    double sigcc(const Event& vertex, const MainEvent&) override { return SigNucElastic(vertex, Z_, A_, M_); }
    double sigccRecon(const Event& recon, const MainEvent&) override { return SigNucElastic(recon, Z_, A_, M_); }

private:
    double Z_{0}, A_{0}, M_{0};
};
}  // namespace

}  // namespace nuclear
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "nuc_elastic", simc::nuclear::NucElasticModel)
