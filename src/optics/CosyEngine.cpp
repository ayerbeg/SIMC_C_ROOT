// src/optics/CosyEngine.cpp — see include/simc/optics/CosyEngine.h.
#include "simc/optics/CosyEngine.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace optics {

// libgcc __powidf2 (what gfortran calls for real*8 ** integer with a non-constant exponent).
double powi(double x, int m) {
    unsigned int n = m < 0 ? -static_cast<unsigned int>(m) : static_cast<unsigned int>(m);
    double y = (n % 2) ? x : 1;
    while (n >>= 1) {
        x = x * x;
        if (n % 2) y = y * x;
    }
    return m < 0 ? 1 / y : y;
}

namespace {

// Fortran: loren.f:loren (l.1-26). Single implementation in setup/KinematicsUtils (P2.2).
void loren(double gam, double bx, double by, double bz, double e, double x, double y, double z, double& ef1,
           double& pxf, double& pyf, double& pzf, double& pf1) {
    const LorenResult r = Loren(gam, bx, by, bz, e, x, y, z);
    ef1 = r.ef1;
    pxf = r.pxf;
    pyf = r.pyf;
    pzf = r.pzf;
    pf1 = r.pf1;
}

// gfortran list-directed integer field (I11 after a one-blank separator).
std::string listInt(long v) {
    char b[32];
    std::snprintf(b, sizeof b, "%11ld", v);
    return b;
}
std::string listReal(double v) {
    char b[48];
    std::snprintf(b, sizeof b, "%25.16f", v);   // approximate (warnings only; not compared in L1)
    return b;
}

// Decay kinematics shared by transp (both halves) and project. Fortran: transp.f:146-184 / 238-274,
// project.f:72-109. m_final_kk: mass given to the 2-pion kaon mode (FORTRAN-BUG D-006: Mk in the first half
// of transp, transp.f:158, kept with fortran_compat; Mpi elsewhere, transp.f:251, project.f:84).
// `rec` (optional) receives the loren inputs and outputs.
struct DecayOut {
    double pxf, pyf, pzf, pf;
};
DecayOut decayProducts(DecayState& d, double m_final_kk, double beta, double gamma, double dxdzs, double dydzs,
                       bool m2BeforeBoost, DecayRecord* rec) {
    using namespace constants;
    IRandom& rng = *d.rng;
    const double rph = rng.uniform() * 2. * pi;
    const double rth1 = rng.uniform() * 2. - 1.;
    const double rth = std::acos(rth1);

    double pr = 0.;
    double m_final = Mmu;   // default
    if (std::fabs(std::sqrt(d.m2) - Mpi) < 2) pr = 29.783;   // pion decay
    if (std::fabs(std::sqrt(d.m2) - Mk) < 2) {                // kaons
        if (rng.uniform() < 0.7) {                            // decay to muon plus neutrino
            pr = 235.5;
        } else {                                              // decay to two pions
            pr = std::sqrt(Mk * Mk / 4. - Mpi * Mpi);
            m_final = m_final_kk;
        }
    }
    if (pr == 0.) {
        // FORTRAN-BUG(D-037 F9): original Fortran shared/project.f:87-91, transp.f:161-165, 254-258: a bare `stop` (exit
        // status 0, no output) for the decay of a particle that is neither a pion nor a kaon. Corrected (default mode)
        // by never checking such a particle for decay (decaysInFlight) and by rejecting doing_decay with a hadron that
        // does not decay at start-up (BuildRunSetup); with fortran_compat the run ends here (simc::FortranStop).
        char b[96];
        std::snprintf(b, sizeof b, "error, cannot decay particle with mass=%8.2f", std::sqrt(d.m2));
        throw FortranStop(b);
    }
    const double er = std::sqrt(m_final * m_final + pr * pr);
    const double pxr = pr * std::sin(rth) * std::cos(rph);
    const double pyr = pr * std::sin(rth) * std::sin(rph);
    const double pzr = pr * std::cos(rth);
    if (m2BeforeBoost) {   // project.f:96-97
        d.m2 = m_final * m_final;
        d.Mh2_final = d.m2;
    }
    const double nrm = std::sqrt(1. + dxdzs * dxdzs + dydzs * dydzs);
    const double bx = -beta * dxdzs / nrm;
    const double by = -beta * dydzs / nrm;
    const double bz = -beta * 1. / nrm;
    DecayOut o{};
    double ef;
    loren(gamma, bx, by, bz, er, pxr, pyr, pzr, ef, o.pxf, o.pyf, o.pzf, o.pf);
    if (!m2BeforeBoost) {  // transp.f:183-184, 273-274
        d.m2 = m_final * m_final;
        d.Mh2_final = d.m2;
    }
    if (rec) {
        rec->gamma = gamma;
        rec->bx = bx;
        rec->by = by;
        rec->bz = bz;
        rec->er = er;
        rec->pxr = pxr;
        rec->pyr = pyr;
        rec->pzr = pzr;
        rec->ef = ef;
        rec->pxf = o.pxf;
        rec->pyf = o.pyf;
        rec->pzf = o.pzf;
        rec->pf = o.pf;
    }
    return o;
}

StepState stepState(const Track& t, const DecayState& d, double pathlen) {
    StepState s;
    s.xs = t.xs;
    s.ys = t.ys;
    s.dxdzs = t.dxdzs;
    s.dydzs = t.dydzs;
    s.dpps = t.dpps;
    s.m2 = d.m2;
    s.ph = d.ph;
    s.pathlen = pathlen;
    s.decdist = d.decdist;
    return s;
}

// Diagnostics of one step (TransportTrace); inactive unless decay->trace is set.
class StepTracer {
public:
    StepTracer(int kind, int cls, double zd, const Track& t, DecayState* d, const double& pathlen, bool check)
        : t_(t), d_(d), pathlen_(pathlen), tr_(d ? d->trace : nullptr) {
        if (!tr_) return;
        ++tr_->n_steps;
        if (check) ++tr_->n_decay_checks;
        if (tr_->steps) {
            rec_.kind = kind;
            rec_.cls = cls;
            rec_.zd = zd;
            rec_.decay_check = check;
            rec_.dflag_in = d->dflag;
            rec_.in = stepState(t, *d, pathlen);
        }
    }
    StepTracer(const StepTracer&) = delete;
    StepTracer& operator=(const StepTracer&) = delete;
    /// The record of a decay in this step (nullptr if not tracing).
    DecayRecord* decay() {
        if (!tr_) return nullptr;
        tr_->decayed = true;
        tr_->dec_check = tr_->n_decay_checks;
        tr_->dec_step = tr_->n_steps;
        decayed_ = true;
        return &tr_->dec;
    }
    ~StepTracer() {
        if (!tr_) return;
        if (decayed_) {
            tr_->dec.dpp_after = t_.dpps;
            tr_->dec.m2_after = d_->m2;
            tr_->dec.p_after = d_->ph;
        }
        if (tr_->steps) {
            rec_.dflag_out = d_->dflag;
            rec_.out = stepState(t_, *d_, pathlen_);
            tr_->steps->push_back(rec_);
        }
    }

private:
    const Track& t_;
    DecayState* d_;
    const double& pathlen_;
    TransportTrace* tr_;
    StepRecord rec_;
    bool decayed_{false};
};

void recordBefore(DecayRecord* r, int kind, int cls, double step_len, const Track& t, const DecayState& d) {
    if (!r) return;
    r->kind = kind;
    r->cls = cls;
    r->step_len = step_len;
    r->dxdz_before = t.dxdzs;
    r->dydz_before = t.dydzs;
    r->dpp_before = t.dpps;
    r->p_before = d.ph;
    r->m2_before = d.m2;
}

void recordAt(DecayRecord* r, const Track& t, const DecayState& d, double pathlen, double z_decay) {
    if (!r) return;
    r->z_decay = z_decay;
    r->x = t.xs;
    r->y = t.ys;
    r->pathlen_at = pathlen;
    r->decdist_at = d.decdist;
}

void recordAfter(DecayRecord* r, const Track& t) {
    if (!r) return;
    r->dxdz_after = t.dxdzs;
    r->dydz_after = t.dydzs;
}

std::string pad132(std::string s) {
    // Fortran read(chan,'(a)') into character*132: truncate / blank-pad to 132.
    if (s.size() > 132) s.resize(132);
    s.resize(132, ' ');
    return s;
}

bool allBlank(const std::string& s) { return s.find_first_not_of(' ') == std::string::npos; }

std::vector<std::string> readLines(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw FortranStop("cannot open COSY file '" + path + "'");
    std::vector<std::string> v;
    std::string l;
    while (std::getline(f, l)) v.push_back(pad132(l));
    return v;
}

// Fortran list-directed read of one real from `s` (read(str_line(10:),*) x).
double listRead(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) throw FortranStop("list-directed read: no value in '" + s + "'");
    std::size_t b = s.find_first_of(" \t,/", a);
    std::string tok = s.substr(a, b == std::string::npos ? std::string::npos : b - a);
    for (auto& c : tok)
        if (c == 'd' || c == 'D' || c == 'q' || c == 'Q') c = 'E';
    return detail::readFortranReal(tok, 0);
}

}  // namespace

namespace detail {

double readFortranReal(const std::string& field, int d) {
    std::string s;
    for (char c : field)
        if (c != ' ' && c != '\t') s += c;   // BLANK='NULL': blanks ignored
    if (s.empty()) return 0.0;
    std::size_t i = 0;
    std::string sign;
    if (s[i] == '+' || s[i] == '-') {
        if (s[i] == '-') sign = "-";
        ++i;
    }
    std::string ip, fp;
    bool point = false;
    for (; i < s.size(); ++i) {
        char c = s[i];
        if (std::isdigit(static_cast<unsigned char>(c))) {
            (point ? fp : ip) += c;
        } else if (c == '.' && !point) {
            point = true;
        } else {
            break;
        }
    }
    long ex = 0;
    if (i < s.size()) {
        char c = s[i];
        if (c == 'E' || c == 'e' || c == 'D' || c == 'd' || c == 'Q' || c == 'q') ++i;
        std::string es = s.substr(i);
        if (es.empty()) throw FortranStop("bad real field '" + field + "'");
        char* end = nullptr;
        ex = std::strtol(es.c_str(), &end, 10);
        if (*end != '\0') throw FortranStop("bad real field '" + field + "'");
    }
    if (!point) ex -= d;   // implied decimal point (Gw.d input editing)
    std::string num = sign + (ip.empty() ? "0" : ip) + "." + fp + "e" + std::to_string(ex);
    return std::strtod(num.c_str(), nullptr);
}

int readFortranInt(const std::string& field) {
    std::string s;
    for (char c : field)
        if (c != ' ') s += c;
    if (s.empty()) return 0;
    char* end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (*end != '\0') throw FortranStop("bad integer field '" + field + "'");
    return static_cast<int>(v);
}

}  // namespace detail

// ---------------------------------------------------------------------------------------------
// project / rotate

bool decaysInFlight(double m2) {
    return std::fabs(std::sqrt(m2) - constants::Mpi) < 2 || std::fabs(std::sqrt(m2) - constants::Mk) < 2;
}

// Fortran: shared/project.f:project (l.1-122)
void project(Track& t, double& x_new, double& y_new, double z_drift, double& pathlen, DecayState* decay) {
    // FORTRAN-BUG(D-037 F9) (see decayProducts): corrected mode checks only pions and kaons for decay.
    const bool check = decay != nullptr && decay->decay_flag && !decay->dflag &&
                       (decay->fortran_compat || decaysInFlight(decay->m2));
    StepTracer tracer(1, 0, z_drift, t, decay, pathlen, check);
    if (!check) {   // particle has already decayed (or no decay), just drift
        pathlen = pathlen + z_drift * std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        x_new = x_new + t.dxdzs * z_drift;
        y_new = y_new + t.dydzs * z_drift;
        return;
    }
    DecayState& d = *decay;
    if (!d.rng) throw std::logic_error("optics::project: decay requested without an IRandom");
    const double p_spec = d.ph / (1. + t.dpps / 100.);
    const double beta = d.ph / std::sqrt(d.ph * d.ph + d.m2);
    const double gamma = 1. / std::sqrt(1. - beta * beta);
    const double dlen = d.ctau * beta * gamma;   // 1/e decay length
    if (z_drift <= 0) std::cout << " drift distance<0:  automatic decay! BAD!\n";
    const double z_decay = -1. * dlen * std::log(1 - d.rng->uniform());

    if (z_decay > z_drift * std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs)) {   // no decay
        d.decdist = d.decdist + z_drift * std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        pathlen = pathlen + z_drift * std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        x_new = x_new + t.dxdzs * z_drift;
        y_new = y_new + t.dydzs * z_drift;
    } else {   // DECAY
        DecayRecord* rec = tracer.decay();
        recordBefore(rec, 1, 0, z_drift, t, d);
        d.dflag = true;
        d.decdist = d.decdist + z_decay;
        pathlen = pathlen + z_decay;
        const double norm_parent = std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        x_new = x_new + t.dxdzs * z_decay / std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        y_new = y_new + t.dydzs * z_decay / std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        recordAt(rec, t, d, pathlen, z_decay);
        DecayOut o = decayProducts(d, constants::Mpi, beta, gamma, t.dxdzs, t.dydzs, true, rec);
        t.dxdzs = o.pxf / o.pzf;
        t.dydzs = o.pyf / o.pzf;
        t.dpps = 100. * (o.pf / p_spec - 1.);
        d.ph = o.pf;
        recordAfter(rec, t);
        // We've already drifted z_decay/sqrt(1+dxdzs**2+dydzs**2) - now do the rest.
        // FORTRAN-BUG(D-037 F7): original Fortran shared/project.f:113:
        // `tmpdrift = z_drift - z_decay/sqrt(1+dxdzs**2+dydzs**2)` evaluated with the slopes AFTER the decay; corrected
        // because the z already covered is z_decay/sqrt(1+s**2) with the PARENT slopes s (project.f:68-69): with the
        // daughter slopes the track ends at z_drift + z_decay*(1/sqrt(1+s_old**2) - 1/sqrt(1+s_new**2)) instead of
        // z_drift (up to 0.28 cm in the P6.F2 captures, 18.9 cm in extreme rows). fortran_compat keeps the Fortran.
        const double tmpdrift = d.fortran_compat
                                    ? z_drift - z_decay / std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs)
                                    : z_drift - z_decay / norm_parent;
        pathlen = pathlen + tmpdrift * std::sqrt(1 + t.dxdzs * t.dxdzs + t.dydzs * t.dydzs);
        x_new = x_new + t.dxdzs * tmpdrift;
        y_new = y_new + t.dydzs * tmpdrift;
    }
}

// FORTRAN-BUG(D-015 F9): original Fortran shared/rotate_haxis.f:41 and rotate_vaxis.f:40:
// `parameter (raddeg=0.017453292)`; corrected because the constant is pi/180 truncated to 8 significant digits
// (relative error 1.1e-8 against pi/180 = 0.017453292519943295); it rotates the HMS dipole and SHMS aperture tests.
// fortran_compat keeps the truncated value.
static double rotateRaddeg(bool fortran_compat) { return fortran_compat ? kRotateRaddeg : constants::pi / 180.; }

// Fortran: shared/rotate_haxis.f:rotate_haxis
void rotateHaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat) {
    const double rotang_rad = rotang * rotateRaddeg(fortran_compat);   // FORTRAN-BUG(D-015 F9) above
    const double tan_th = std::tan(rotang_rad);
    const double sin_th = std::sin(rotang_rad);
    const double cos_th = std::cos(rotang_rad);
    const double alpha = t.dxdzs;
    const double beta = t.dydzs;
    const double alpha_p = (alpha + tan_th) / (1. - alpha * tan_th);
    const double beta_p = beta / (cos_th - alpha * sin_th);
    const double xi = xp0;
    xp0 = xi * (cos_th + alpha_p * sin_th);
    yp0 = yp0 + xi * beta_p * sin_th;
}

// Fortran: shared/rotate_vaxis.f:rotate_vaxis
void rotateVaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat) {
    const double rotang_rad = rotang * rotateRaddeg(fortran_compat);   // FORTRAN-BUG(D-015 F9) above
    const double tan_th = std::tan(rotang_rad);
    const double sin_th = std::sin(rotang_rad);
    const double cos_th = std::cos(rotang_rad);
    const double alpha = t.dydzs;
    const double beta = t.dxdzs;
    const double alpha_p = (alpha + tan_th) / (1. - alpha * tan_th);
    const double beta_p = beta / (cos_th - alpha * sin_th);
    const double yi = yp0;
    yp0 = yi * (cos_th + alpha_p * sin_th);
    xp0 = xp0 + yi * beta_p * sin_th;
}

// ---------------------------------------------------------------------------------------------
// CosyEngine

CosyEngine::CosyEngine() : log_(&std::cout) {}

// Fortran: transp.f:300-314
std::string CosyEngine::defaultForwardFile(int spectr) {
    switch (spectr) {
        case 1: return "hms/forward_cosy.dat";
        case 2: return "sos/forward_cosy.dat";
        case 3: return "hrsr/hrs_forward_cosy.dat";
        case 4: return "hrsl/hrs_forward_cosy.dat";
        case 5: return "shms/shms_forward.dat";
        case 6: throw FortranStop("LSA tune for SHMS no longer used! (transp.f:312)");
        default: throw FortranStop("transp_init: unknown spectrometer " + std::to_string(spectr));
    }
}

// Fortran: mc_hms_recon.f:73, mc_sos_recon.f:73, mc_shms_recon.f:79
std::string CosyEngine::defaultReconFile(int spectr) {
    if (spectr == 1) return "hms/recon_cosy.dat";
    if (spectr == 2) return "sos/recon_cosy.dat";
    if (spectr == 5) return "shms/shms_recon.dat";
    // (mc_shms_recon.f also STOPs for the retired SHMS LSA tune; the input rejects that arm code, D-032.)
    throw FortranStop("recon for spectrometer " + std::to_string(spectr) + " not ported (HMS 1, SOS 2, SHMS 5 only)");
}

int CosyEngine::initForward(int spectr, const std::string& simcDir) {
    return loadForward(spectr, simcDir + "/" + defaultForwardFile(spectr));
}
int CosyEngine::initRecon(int spectr, const std::string& simcDir) {
    return loadRecon(spectr, simcDir + "/" + defaultReconFile(spectr));
}

// Fortran: transp.f:transp_init / transp_init_file (l.285-478)
int CosyEngine::loadForward(int spectr, const std::string& path) {
    const std::vector<std::string> lines = readLines(path);
    std::size_t idx = 0;
    auto next = [&](std::string& s) {
        if (idx >= lines.size()) return false;
        s = lines[idx++];
        return true;
    };
    // The Fortran `length` array is static and only set from !LENGTH: comments; a fresh load starts at 0.
    std::vector<double> length(kMaxClass + 2, 0.0);
    std::vector<ForwardClass> cls;

    // Strip away header.
    std::string str = pad132("!");
    while (str[0] == '!') {
        if (!next(str)) throw FortranStop("TRANSP_INIT: end of file in header of " + path);
        if (str.compare(0, 8, "!LENGTH:") == 0) length[1] = 100. * listRead(str.substr(9));
    }

    int n_classes = 0;
    for (;;) {
        const int kk = n_classes + 1;
        if (kk > kMaxClass) throw FortranStop("TRANSP_INIT: too many transformations!");
        ForwardClass c;   // adrift = .true., driftdist = 0
        while (str.compare(0, 4, " ---") != 0) {
            if (static_cast<int>(c.terms.size()) + 1 > kMaxElements) throw FortranStop("TRANSP_INIT: too many COSY terms!");
            // format(1x,5g14.7,1x,6i1): coeff X A Y B dL; exponents x a y b TOF delta
            ForwardTerm t;
            for (int i = 0; i < 5; ++i) t.coeff[i] = detail::readFortranReal(str.substr(1 + 14 * i, 14), 7);
            int e[6];
            for (int j = 0; j < 6; ++j) e[j] = detail::readFortranInt(str.substr(72 + j, 1));
            t.expon = {e[0], e[1], e[2], e[3], e[5]};
            const int idummy = e[4];
            if (idummy != 0) {   // Ignore time-of-flight term.
                if (t.coeff[0] != 0 || t.coeff[1] != 0 || t.coeff[2] != 0 || t.coeff[3] != 0)
                    throw FortranStop("TRANSP_INIT: non-zero TOF terms!");
            } else {
                c.terms.push_back(t);
                // Check to see if element is inconsistent with a field-free drift.
                if (c.adrift) {
                    const int e1 = t.expon[0], e2 = t.expon[1], e3 = t.expon[2], e4 = t.expon[3];
                    const double c1 = t.coeff[0], c2 = t.coeff[1], c3 = t.coeff[2], c4 = t.coeff[3];
                    const double csum = std::fabs(c1) + std::fabs(c2) + std::fabs(c3) + std::fabs(c4);
                    const int order = e1 + e2 + e3 + e4;
                    if (order == 1) {
                        if (e1 == 1) {
                            if (std::fabs(c1 - 1.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c2 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c3 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c4 - 0.e0) > kCoeffMin) c.adrift = false;
                        } else if (e2 == 1) {
                            c.driftdist = 1000.e0 * c1;   // drift distance in cm.
                            if (std::fabs(c2 - 1.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c3 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c4 - 0.e0) > kCoeffMin) c.adrift = false;
                        } else if (e3 == 1) {
                            if (std::fabs(c1 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c2 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c3 - 1.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c4 - 0.e0) > kCoeffMin) c.adrift = false;
                        } else if (e4 == 1) {
                            if (std::fabs(c1 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c2 - 0.e0) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c.driftdist - 1000.e0 * c3) > kCoeffMin) c.adrift = false;
                            if (std::fabs(c4 - 1.e0) > kCoeffMin) c.adrift = false;
                        }
                    } else {   // if order.ne.1
                        if (std::fabs(csum) > kCoeffMin) c.adrift = false;
                    }
                }
            }
            if (!next(str)) throw FortranStop("TRANSP_INIT: end of file inside a class of " + path);
        }
        n_classes = kk;
        c.length = length[kk];
        if (c.adrift && std::fabs(c.driftdist - c.length) > 0.01 && log_) {
            *log_ << " PROBLEM WITH TRANSFORMATION #" << listInt(kk) << " for spectrometer #" << listInt(spectr) << "\n";
            *log_ << " Appears to be a pure drift, driftdist=" << listReal(c.driftdist) << "\n";
            *log_ << " But the comments say that the length =" << listReal(c.length) << "\n";
        }
        cls.push_back(c);
        // Read lines until a non-blank, non-comment non-terminal line is found.
        bool eof = false;
        for (;;) {
            if (!next(str)) {
                eof = true;
                break;
            }
            if (str.compare(0, 8, "!LENGTH:") == 0) length[kk + 1] = 100. * listRead(str.substr(9));
            if (str[0] == '!' || str.compare(0, 4, " ---") == 0 || allBlank(str)) continue;
            break;
        }
        if (eof) break;
    }
    fwd_[spectr] = std::move(cls);
    return n_classes;
}

// Fortran: mc_hms_recon.f:64-93 / mc_shms_recon.f:66-119 (firsttime block)
int CosyEngine::loadRecon(int spectr, const std::string& path) {
    const std::vector<std::string> lines = readLines(path);
    std::size_t idx = 0;
    auto next = [&](std::string& s) {
        if (idx >= lines.size()) return false;
        s = lines[idx++];
        return true;
    };
    std::string line = pad132("!");
    while (line[0] == '!') {
        if (!next(line)) throw FortranStop("recon: end of file in header of " + path);
    }
    std::vector<ReconTerm> terms;
    int max_order = 0;
    while (line.compare(0, 4, " ---") != 0) {
        if (static_cast<int>(terms.size()) + 1 > kMaxElements) throw FortranStop("WCRECON: too many COSY terms!");
        // format(1x,4g16.9,1x,5i1)
        ReconTerm t;
        for (int i = 0; i < 4; ++i) t.coeff[i] = detail::readFortranReal(line.substr(1 + 16 * i, 16), 9);
        for (int j = 0; j < 5; ++j) t.expon[j] = detail::readFortranInt(line.substr(66 + j, 1));
        terms.push_back(t);
        if (!next(line)) throw FortranStop("recon: end of file before ' ---' in " + path);
        max_order = std::max(max_order, t.expon[0] + t.expon[1] + t.expon[2] + t.expon[3] + t.expon[4]);
    }
    rec_[spectr] = std::move(terms);
    recMaxOrder_[spectr] = max_order;
    return static_cast<int>(rec_[spectr].size());
}

int CosyEngine::nClasses(int spectr) const {
    auto it = fwd_.find(spectr);
    return it == fwd_.end() ? 0 : static_cast<int>(it->second.size());
}

const CosyEngine::ForwardClass& CosyEngine::forwardClass(int spectr, int cls) const {
    auto it = fwd_.find(spectr);
    if (it == fwd_.end()) throw std::out_of_range("CosyEngine: no forward map for spectrometer " + std::to_string(spectr));
    if (cls < 1 || cls > static_cast<int>(it->second.size()))
        throw std::out_of_range("CosyEngine: class " + std::to_string(cls) + " out of range for spectrometer " +
                                std::to_string(spectr));
    return it->second[static_cast<std::size_t>(cls - 1)];
}

const std::vector<CosyEngine::ReconTerm>& CosyEngine::reconTerms(int spectr) const {
    auto it = rec_.find(spectr);
    if (it == rec_.end()) throw std::out_of_range("CosyEngine: no recon map for spectrometer " + std::to_string(spectr));
    return it->second;
}

int CosyEngine::reconMaxOrder(int spectr) const {
    auto it = recMaxOrder_.find(spectr);
    return it == recMaxOrder_.end() ? 0 : it->second;
}

// Fortran: shared/transp.f:transp (l.1-281)
double CosyEngine::transp(int spectr, int cls, Track& t, double zd, double& pathlen, DecayState* decay) {
    const ForwardClass& c = forwardClass(spectr, cls);
    if (c.length > 0.01 && std::fabs(c.length - zd) > 0.01 && log_) {
        *log_ << " PROBLEM WITH TRANSFORMATION #" << listInt(cls) << " for spectrometer #" << listInt(spectr) << "\n";
        *log_ << " Central path length passed to transp =" << listReal(zd) << "\n";
        *log_ << " But the comments say that the length =" << listReal(c.length) << "\n";
    }

    // Check for decay (first half of the element).
    // FORTRAN-BUG(D-037 F9) (see decayProducts): corrected mode checks only pions and kaons for decay.
    const bool check = decay != nullptr && decay->decay_flag && !decay->dflag &&
                       (fortranCompat_ || decaysInFlight(decay->m2));
    StepTracer tracer(2, cls, zd, t, decay, pathlen, check);
    DecayRecord* rec = nullptr;
    double p_spec = 0, beta = 0, gamma = 0, z_decay = 0;
    if (check) {
        DecayState& d = *decay;
        if (!d.rng) throw std::logic_error("optics::transp: decay requested without an IRandom");
        p_spec = d.ph / (1. + t.dpps / 100.);
        beta = d.ph / std::sqrt(d.ph * d.ph + d.m2);
        gamma = 1. / std::sqrt(1. - beta * beta);
        const double dlen = d.ctau * beta * gamma;
        if (zd <= 0 && log_) *log_ << " transp distance<0:  automatic decay! BAD!\n";
        z_decay = -1. * dlen * std::log(1 - d.rng->uniform());
        if (z_decay <= zd / 2) {   // decay in first half
            // FORTRAN-DIFF-CANDIDATE(D-037 F8): the first-half decay is applied at the element ENTRANCE (position
            // unchanged, then the COSY map transports the daughter) and the second-half decay at the element EXIT
            // with the parent's beta/gamma; decdist gets z_decay while pathlen always gets zd+delta_z. Question: is
            // this approximation intended (the map cannot be split)? Kept in both modes.
            rec = tracer.decay();
            recordBefore(rec, 2, cls, zd, t, d);
            d.dflag = true;
            d.decdist = d.decdist + z_decay;
            recordAt(rec, t, d, pathlen, z_decay);
            // FORTRAN-BUG(D-006): original Fortran shared/transp.f:158: `m_final = Mk` for the K -> pi pi mode of a decay
            // in the first half of an element; corrected because the decay product is a pion (pr uses Mpi, and the
            // second half transp.f:251 and project.f:84 set m_final = Mpi). fortran_compat keeps Mk.
            DecayOut o = decayProducts(d, fortranCompat_ ? constants::Mk : constants::Mpi, beta, gamma, t.dxdzs, t.dydzs,
                                       false, rec);
            t.dxdzs = o.pxf / o.pzf;
            t.dydzs = o.pyf / o.pzf;
            t.dpps = 100. * (o.pf / p_spec - 1.);
            d.ph = o.pf;
            recordAfter(rec, t);
        }
    }

    // Pack local copy of input coordinates.
    const double ray[5] = {t.xs, t.dxdzs * 1000.e0, t.ys, t.dydzs * 1000.e0, t.dpps};
    double sum[5] = {0., 0., 0., 0., 0.};
    for (const ForwardTerm& e : c.terms) {
        double term = 1.0e0;
        for (int j = 0; j < 5; ++j) {
            double temp = 1.0e0;
            if (e.expon[j] != 0) temp = powi(ray[j], e.expon[j]);
            term = term * temp;
        }
        sum[0] = sum[0] + term * e.coeff[0];   // NEW "X"
        sum[1] = sum[1] + term * e.coeff[1];   // NEW "A"
        sum[2] = sum[2] + term * e.coeff[2];   // NEW "Y"
        sum[3] = sum[3] + term * e.coeff[3];   // NEW "B"
        sum[4] = sum[4] + term * e.coeff[4];   // NEW "dL"
    }
    t.xs = sum[0];
    t.dxdzs = sum[1] / 1000.e0;
    t.ys = sum[2];
    t.dydzs = sum[3] / 1000.e0;
    const double delta_z = -sum[4];

    // Check for decay in 2nd half of element, which is applied AFTER transporting.
    if (check && !decay->dflag) {
        DecayState& d = *decay;
        if (z_decay > zd + delta_z) {   // don't decay within magnet.
            d.decdist = d.decdist + (zd + delta_z);
        } else {                        // apply decay at end/magnet.
            rec = tracer.decay();
            recordBefore(rec, 3, cls, zd, t, d);
            d.dflag = true;
            d.decdist = d.decdist + z_decay;
            recordAt(rec, t, d, pathlen, z_decay);
            DecayOut o = decayProducts(d, constants::Mpi, beta, gamma, t.dxdzs, t.dydzs, false, rec);
            t.dxdzs = o.pxf / o.pzf;
            t.dydzs = o.pyf / o.pzf;
            t.dpps = 100. * (o.pf / p_spec - 1.);
            d.ph = o.pf;
            recordAfter(rec, t);
        }
    }

    pathlen = pathlen + (zd + delta_z);
    return delta_z;
}

// Fortran: hms/mc_hms_recon.f (l.95-129), shms/mc_shms_recon.f (l.121-155), sos/mc_sos_recon.f (l.104-141) — identical
// arithmetic, except the zero guard: HMS/SHMS replace only hut(5) (fry) by 1e-30, SOS all five hut values
// (mc_sos_recon.f:120-122).
void CosyEngine::recon(int spectr, const Track& t, double fry, double& delta_p, double& delta_t, double& delta_phi,
                       double& y_tgt) const {
    const std::vector<ReconTerm>& terms = reconTerms(spectr);
    double sum[4] = {0., 0., 0., 0.};
    double hut[5];
    hut[0] = t.xs / 100.;   // cm --> m
    hut[1] = t.dxdzs;       // slope
    hut[2] = t.ys / 100.;   // cm --> m
    hut[3] = t.dydzs;       // slope
    hut[4] = fry / 100.;    // vert. position at target (cm --> m)
    if (spectr == 2) {
        // mc_sos_recon.f:120-122: do i=1,5; if (abs(hut(i)).le.1.e-30) hut(i)=1.e-30
        for (int i = 0; i < 5; ++i)
            if (std::fabs(hut[i]) <= 1.e-30) hut[i] = 1.e-30;
    } else {
        if (std::fabs(hut[4]) <= 1.e-30) hut[4] = 1.e-30;   // FORTRAN-QUIRK (D-015): avoid 0.0**0
    }
    for (const ReconTerm& e : terms) {
        const double term = powi(hut[0], e.expon[0]) * powi(hut[1], e.expon[1]) * powi(hut[2], e.expon[2]) *
                            powi(hut[3], e.expon[3]) * powi(hut[4], e.expon[4]);
        sum[0] = sum[0] + term * e.coeff[0];
        sum[1] = sum[1] + term * e.coeff[1];
        sum[2] = sum[2] + term * e.coeff[2];
        sum[3] = sum[3] + term * e.coeff[3];
    }
    delta_phi = sum[0];        // slope
    y_tgt = sum[1] * 100.;     // m --> cm
    delta_t = sum[2];          // slope
    delta_p = sum[3] * 100.;   // percent deviation
}

}  // namespace optics
}  // namespace simc
