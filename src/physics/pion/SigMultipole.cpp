// src/physics/pion/SigMultipole.cpp
//
// Fortran: physics_pion.f:sig_multipole (l.197-398). Not called by SIMC (its call in peepi is commented
// out, physics_pion.f:157-163); ported for the L1 table sig_multipole.
//
// FORTRAN-DIFF-CANDIDATE(D-038 F12): COMPLEX*8 amplitudes under -fdefault-real-8 (kept in both modes; Question: is
// single precision intended?). The Fortran kind rules are applied operation
// by operation:
//   * complex(4) op complex(4)           -> computed in single precision (float);
//   * real(8) * complex(4), complex(4)/real(8), complex(4) + complex(8) ... -> COMPLEX(8) (the complex(4) operand is
//     widened exactly), and the assignment to a COMPLEX*8 variable rounds each part to float;
//   * a real operand is not turned into a full complex product: gfortran -O lowers real*complex and complex/real
//     component-wise (GCC tree-complex "only real" lattice), complex*complex with the plain formula (-fcx-fortran-rules:
//     no NaN recovery): (a+ib)(c+id) = (ac-bd) + i(ad+bc);
//   * REAL(complex(4)) is real(4); CMPLX(x,y) without KIND is the default COMPLEX(8) (rounded to COMPLEX*8 on
//     assignment).
// std::complex<float> is only the storage; the arithmetic is written out with the helpers below so that no library
// complex operator (C99 Annex G NaN handling, __mulsc3) is involved.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "simc/physics/pion/PionCrossSections.h"

namespace simc {
namespace pion {

namespace {
using cf = std::complex<float>;
using cd = std::complex<double>;

inline cd W(cf a) { return cd(static_cast<double>(a.real()), static_cast<double>(a.imag())); }   // widen (exact)
inline cf N(cd a) { return cf(static_cast<float>(a.real()), static_cast<float>(a.imag())); }     // assign to COMPLEX*8
inline cd mulR(double k, cd a) { return cd(k * a.real(), k * a.imag()); }                       // real(8)*complex(8)
inline cd divR(cd a, double k) { return cd(a.real() / k, a.imag() / k); }                       // complex(8)/real(8)
inline cd add(cd a, cd b) { return cd(a.real() + b.real(), a.imag() + b.imag()); }
inline cd sub(cd a, cd b) { return cd(a.real() - b.real(), a.imag() - b.imag()); }
inline cd neg(cd a) { return cd(-a.real(), -a.imag()); }
inline cf addF(cf a, cf b) { return cf(a.real() + b.real(), a.imag() + b.imag()); }
inline cf subF(cf a, cf b) { return cf(a.real() - b.real(), a.imag() - b.imag()); }
inline cf mulF(cf a, cf b) {
    return cf(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}
inline cf conjF(cf a) { return cf(a.real(), -a.imag()); }

/// Fortran F edit descriptor on one record (blanks ignored, PAD with blanks, implied decimals).
double fField(const std::string& rec, std::size_t pos, std::size_t w, int d) {
    std::string s;
    for (std::size_t i = pos; i < pos + w && i < rec.size(); ++i)
        if (rec[i] != ' ') s += rec[i];
    if (s.empty()) return 0.0;
    double v = std::strtod(s.c_str(), nullptr);
    if (s.find('.') == std::string::npos && s.find_first_of("eEdD") == std::string::npos) v = v / std::pow(10.0, d);
    return v;
}
}  // namespace

// Fortran: physics_pion.f:255-286 (OPEN(3,FILE='MULTIPOL.file'), formats 202-204 l.394-396)
bool SigMultipole::load(const std::string& path, std::string* err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        if (err) *err = "sig_multipole: cannot open " + path;
        return false;
    }
    std::vector<std::string> recs;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        recs.push_back(line);
    }
    std::size_t r = 0;   // next record
    auto next = [&](std::string& out) {
        if (r >= recs.size()) return false;
        out = recs[r++];
        return true;
    };
    std::string rec;
    for (int q2c = 1; q2c <= 8; ++q2c) {
        // 202 format(/11X,f5.1/): skip a record, read qfm, skip a record
        std::string dummy;
        if (!next(dummy) || !next(rec) || !next(dummy)) {
            if (err) *err = "sig_multipole: " + path + " too short";
            return false;
        }
        const double qfm = fField(rec, 11, 5, 1);   // q**2 in fm**-2
        q2_table_[q2c - 1] = qfm * 0.0389379;        // -> GeV**2
        for (int nuc = 1; nuc <= 10; ++nuc) {
            // 203 format(11X,f5.0)
            if (!next(rec)) {
                if (err) *err = "sig_multipole: " + path + " too short";
                return false;
            }
            nu_table_[nuc - 1] = fField(rec, 11, 5, 0);   // overwritten by every Q2 block (the last one wins)
            // 204 format(6(/9X,7f8.3)): the read starts on the next record (the column header "E0+ M1- ..."), which
            // the leading '/' leaves unread; each group then reads one record (RE/IM 3/2, 1/2, 0).
            double atemp[42];
            if (!next(dummy)) {
                if (err) *err = "sig_multipole: " + path + " too short";
                return false;
            }
            for (int g = 0; g < 6; ++g) {
                if (!next(rec)) {
                    if (err) *err = "sig_multipole: " + path + " too short";
                    return false;
                }
                for (int k = 0; k < 7; ++k) atemp[g * 7 + k] = fField(rec, 9 + 8 * static_cast<std::size_t>(k), 8, 3);
            }
            for (int m = 1; m <= 7; ++m) {
                const int it1 = m + 7, it2 = m + 14, it3 = m + 21, it4 = m + 28, it5 = m + 35;
                // AT(m,3,..) = CMPLX(ATEMP(m), ATEMP(IT1)): default (8-byte-real) complex, rounded to COMPLEX*8
                at_[m - 1][2][q2c - 1][nuc - 1] = N(cd(atemp[m - 1], atemp[it1 - 1]));
                at_[m - 1][1][q2c - 1][nuc - 1] = N(cd(atemp[it2 - 1], atemp[it3 - 1]));
                at_[m - 1][0][q2c - 1][nuc - 1] = N(cd(atemp[it4 - 1], atemp[it5 - 1]));
            }
        }
    }
    loaded_ = true;
    return true;
}

// Fortran: physics_pion.f:sig_multipole (l.197-398)
double SigMultipole::operator()(double& nu, double efer, double qstar, double tfcos, double ppicm, double s_gev,
                                double thetacm, double phicm, double eps, int final_state, double q2_gev,
                                double mtar_gev, double pfer, double mh, std::ostream* log) {
    if (nu >= 500.0) {
        if (log) *log << " LOOK OUT CHEESY POOF BREATH!\n W IS TOO HIGH FOR THIS MODEL\n SETTING NU to 499.9\n";
        nu = 499.9;
    }

    const double sintcm = std::sin(thetacm);
    const double costcm = std::cos(thetacm);
    const double sinto2cm = std::sin(thetacm / 2.);
    const double costo2cm = std::cos(thetacm / 2.);

    // l.248-259: first call (the file itself is read by load())
    if (first_) {
        first_ = false;
        root_two_ = std::sqrt(2.);
        mpi_to_microbarn_ = 197.327 / mh;   // SAVEd: the mh of the first call
        cleb_[0][0] = root_two_;
        cleb_[0][1] = root_two_ / 3.;
        cleb_[0][2] = -root_two_ / 3.;
        cleb_[1][0] = root_two_;
        cleb_[1][1] = -root_two_ / 3.;
        cleb_[1][2] = root_two_ / 3.;
    }
    auto CLEB = [&](int fs, int iso) { return cleb_[fs - 1][iso - 1]; };

    // l.289-321: interpolate the multipole amplitudes (strictly inside a cell; otherwise the SAVEd amplitudes stay)
    double Q2_high, Q2_low, nu_high, nu_low, deltaQ2, delta_nu;
    for (int q2c = 1; q2c <= 7; ++q2c) {
        Q2_high = 0.;
        Q2_low = 0.;
        nu_high = 0.;
        nu_low = 0.;
        if (q2_gev > q2_table_[q2c - 1] && q2_gev < q2_table_[q2c]) {
            Q2_high = q2_table_[q2c];
            Q2_low = q2_table_[q2c - 1];
            deltaQ2 = (Q2_high - Q2_low);
            for (int nuc = 1; nuc <= 9; ++nuc) {
                if (nu > nu_table_[nuc - 1] && nu < nu_table_[nuc]) {
                    nu_high = nu_table_[nuc];
                    nu_low = nu_table_[nuc - 1];
                    delta_nu = nu_high - nu_low;
                    for (int m = 1; m <= 7; ++m) {
                        amplitude_[m - 1] = cf(0.f, 0.f);
                        for (int iso = 1; iso <= 3; ++iso) {
                            const cf A1 = at_[m - 1][iso - 1][q2c - 1][nuc - 1];
                            const cf A2 = at_[m - 1][iso - 1][q2c][nuc - 1];
                            const cf A3 = at_[m - 1][iso - 1][q2c - 1][nuc];
                            const cf A4 = at_[m - 1][iso - 1][q2c][nuc];
                            // A12 = (A1*(Q2_high-q2_gev) + A2*(q2_gev-Q2_low))/deltaQ2
                            const cf A12 =
                                N(divR(add(mulR(Q2_high - q2_gev, W(A1)), mulR(q2_gev - Q2_low, W(A2))), deltaQ2));
                            const cf A34 =
                                N(divR(add(mulR(Q2_high - q2_gev, W(A3)), mulR(q2_gev - Q2_low, W(A4))), deltaQ2));
                            // A1234 = (A12*(nu_high-nu)+A34*(nu-nu_low))/delta_nu
                            const cf A1234 =
                                N(divR(add(mulR(nu_high - nu, W(A12)), mulR(nu - nu_low, W(A34))), delta_nu));
                            // amplitude = amplitude + CLEB(final_state,isospin)*A1234
                            amplitude_[m - 1] = N(add(W(amplitude_[m - 1]), mulR(CLEB(final_state, iso), W(A1234))));
                        }
                        // convert multipoles to microb^.5
                        amplitude_[m - 1] = N(mulR(mpi_to_microbarn_, W(amplitude_[m - 1])));
                    }
                }
            }
        }
    }
    (void)nu_high;
    (void)nu_low;
    const cf* a = amplitude_;   // a[k-1] = amplitude(k)

    // l.325-334: CGLN amplitudes
    const cf gf1 = N(add(W(a[0]), mulR(3. * costcm, W(addF(a[3], a[2])))));
    const cf gf2 = N(add(mulR(2., W(a[3])), W(a[1])));
    const cf gf3 = N(mulR(3., W(subF(a[2], a[3]))));
    const cf gf4(0.f, 0.f);
    const cf gf7 = N(sub(W(a[5]), mulR(2., W(a[6]))));
    const cf gf8 = N(add(W(a[4]), mulR(6. * costcm, W(a[6]))));

    // l.337-342: Walker helicity amplitudes
    const cf H5 = N(neg(mulR(costo2cm, W(addF(gf7, gf8)))));
    const cf H6 = N(neg(mulR(sinto2cm, W(subF(gf7, gf8)))));
    const cf H4 = N(divR(add(mulR(2. * sinto2cm, W(addF(gf1, gf2))), mulR(costo2cm * sintcm, W(addF(gf3, gf4)))),
                         root_two_));
    const cf H2 = N(divR(add(neg(mulR(2. * costo2cm, W(subF(gf1, gf2)))), mulR(sinto2cm * sintcm, W(subF(gf3, gf4)))),
                         root_two_));
    const cf H1 = N(divR(neg(mulR(costo2cm * sintcm, W(addF(gf3, gf4)))), root_two_));
    const cf H3 = N(divR(neg(mulR(sinto2cm * sintcm, W(subF(gf3, gf4)))), root_two_));

    // l.353-358: Bartl helicity amplitudes
    const cf hnperp = N(divR(W(addF(H4, H1)), root_two_));
    const cf hfpar = N(divR(W(addF(H3, H2)), root_two_));
    const cf hnlong = H5;
    const cf hflong = H6;
    const cf hnpar = N(divR(W(subF(H4, H1)), root_two_));
    const cf hfperp = N(divR(W(subF(H3, H2)), root_two_));

    // l.360-366: all in single precision
    const cf XT = addF(addF(addF(mulF(hnperp, conjF(hnperp)), mulF(hfpar, conjF(hfpar))), mulF(hnpar, conjF(hnpar))),
                       mulF(hfperp, conjF(hfperp)));
    const cf XL = addF(mulF(hnlong, conjF(hnlong)), mulF(hflong, conjF(hflong)));
    const cf XTT = subF(subF(addF(mulF(hnpar, conjF(hnpar)), mulF(hfperp, conjF(hfperp))), mulF(hnperp, conjF(hnperp))),
                        mulF(hfpar, conjF(hfpar)));
    const cf XLT = addF(mulF(hnlong, conjF(hnpar)), mulF(hflong, conjF(hfpar)));

    const double F = 2. * ((efer - pfer * tfcos) / 1000.) * std::sqrt(s_gev) * (ppicm / 1000.) /
                     (s_gev - mtar_gev * mtar_gev) / mtar_gev;

    // sigt = F/2.*XT: COMPLEX(8) product, the assignment to a real takes the real part
    const double sigt = (F / 2.) * static_cast<double>(XT.real());
    const double sigl = F * static_cast<double>(XL.real());
    const double sigtt = (F / 2.) * static_cast<double>(XTT.real());
    double siglt = 2. * F * static_cast<double>(XLT.real());   // REAL(XLT) is real(4)
    siglt = siglt / 2.;   // Henk's convention (DJG)

    const double sig219 = (sigt + eps * sigl + eps * std::cos(2. * phicm) * sigtt +
                           std::sqrt(2.0 * eps * (1. + eps)) * std::cos(phicm) * siglt) /
                          1.e0;

    // dsig/dOmega_cm -> dsig/dt/dphi_cm with dt/dcos_cm = 2*ppicm*qcm
    return sig219 / ppicm / qstar / 2.;
}

}  // namespace pion
}  // namespace simc
