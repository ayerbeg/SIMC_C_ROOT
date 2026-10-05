// src/physics/nuclear/SpectralFunctionTables.cpp — see include/simc/physics/nuclear/SpectralFunctionTables.h.
// Fortran: simc_gfortran/sf_lookup.f (sf_lookup_init, sf_lookup_diff, sf_lookup, generate_em) and event.f:complete_main
// l.1402-1430. Operation order as in the Fortran (-fdefault-real-8: every literal is a double; x**2 -> x*x).
#include "simc/physics/nuclear/SpectralFunctionTables.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace nuclear {

namespace {

/// Just enough of the Fortran list-directed `read(u,*)` for these files: every read starts on a new record, takes
/// values across records if needed and discards the rest of the last record it touched.
class RecordReader {
public:
    explicit RecordReader(std::istream& in) : in_(in) {}
    bool read(std::vector<double>& v, std::size_t n) {
        v.clear();
        std::string line;
        while (v.size() < n) {
            if (!std::getline(in_, line)) return false;
            for (char& c : line)
                if (c == ',') c = ' ';
            std::istringstream is(line);
            std::string tok;
            while (v.size() < n && is >> tok) {
                for (char& c : tok)
                    if (c == 'd' || c == 'D') c = 'E';
                char* end = nullptr;
                const double x = std::strtod(tok.c_str(), &end);
                if (end == tok.c_str() || *end != '\0') return false;
                v.push_back(x);
            }
        }
        return true;
    }

private:
    std::istream& in_;
};

std::string listInt(long i) {   // gfortran list-directed integer*4: I12
    char b[32];
    std::snprintf(b, sizeof b, "%12ld", i);
    return b;
}

}  // namespace

// =====================================================================================================================
// Fortran: sf_lookup.f:sf_lookup_init (l.1-81)
// =====================================================================================================================
bool SfTable::load(const std::string& path, bool protonflag, std::vector<std::string>* msgs, std::string* err) {
    std::ifstream f(path);
    if (!f) {
        if (err) *err = "sf_lookup_init: cannot open " + path;
        return false;
    }
    RecordReader rd(f);
    std::vector<double> v;
    if (!rd.read(v, 2)) {   // read (55,*) numPm , numEm
        if (err) *err = "sf_lookup_init: cannot read numPm, numEm from " + path;
        return false;
    }
    const int nP = static_cast<int>(v[0]), nE = static_cast<int>(v[1]);
    if (nP < 2 || nE < 2 || nP > nPmmax || nE > nEmmax) {
        if (err) *err = "sf_lookup_init: numPm/numEm out of range (sf_lookup.inc nPmmax=100, nEmmax=200) in " + path;
        return false;
    }
    file_ = path;
    numPm_ = nP;
    numEm_ = nE;
    Pmval_.assign(nPmmax + 1, 0.0);
    dPm_.assign(nPmmax + 1, 0.0);
    Emval_.assign(nEmmax + 1, 0.0);
    dEm_.assign(nEmmax + 1, 0.0);
    sfnorm_.assign(nPmmax + 1, 0.0);
    sfval_.assign(static_cast<std::size_t>((nPmmax + 1) * (nEmmax + 1)), 0.0);
    for (int iPm = 1; iPm <= numPm_; ++iPm) {
        for (int iEm = 1; iEm <= numEm_; ++iEm) {
            if (!rd.read(v, 6)) {   // read (55,*) tmpPm,tmpEm,tmpsfp,tmpsfn,tmpdPm,tmpdEm
                if (err) *err = "sf_lookup_init: read error in " + path;
                numPm_ = numEm_ = 0;
                return false;
            }
            const double tmpPm = v[0], tmpEm = v[1], tmpsfp = v[2], tmpsfn = v[3], tmpdPm = v[4], tmpdEm = v[5];
            sfval_[idx(iEm, iPm)] = protonflag ? tmpsfp : tmpsfn;
            if (iEm == 1) {
                Pmval_[iPm] = tmpPm;
                dPm_[iPm] = tmpdPm;
            } else if (std::abs(tmpPm - Pmval_[iPm]) > 0.01) {
                if (msgs) {
                    msgs->push_back(" Bad Pm value, line number is approx" + listInt(numEm_ * (iPm - 1) + iEm));
                    msgs->push_back(" iPm,oldPm,newPm=" + listInt(iPm) + fortran::listReal(Pmval_[iPm]) +
                                    fortran::listReal(tmpPm));
                }
            }
            if (iPm == 1) {
                Emval_[iEm] = tmpEm;
                dEm_[iEm] = tmpdEm;
            } else if (std::abs(tmpEm - Emval_[iEm]) > 0.01) {
                if (msgs) {
                    msgs->push_back(" Bad Em value, line number is approx" + listInt(numEm_ * (iPm - 1) + iEm));
                    msgs->push_back(" iEm,oldEm,newEm=" + listInt(iEm) + fortran::listReal(Emval_[iEm]) +
                                    fortran::listReal(tmpEm));
                }
            }
        }
    }
    // l.63-78
    sftotnorm_ = 0.0;
    for (int iPm = 1; iPm <= numPm_; ++iPm) {
        sfnorm_[iPm] = 0.0;
        for (int iEm = 1; iEm <= numEm_; ++iEm) {
            sfnorm_[iPm] = sfnorm_[iPm] + sfval_[idx(iEm, iPm)];
            sftotnorm_ = sftotnorm_ + sfval_[idx(iEm, iPm)];
        }
    }
    if (msgs) msgs->push_back(" Uncorrected S.F. has normalization= " + fortran::listReal(sftotnorm_));
    for (int iPm = 1; iPm <= numPm_; ++iPm)
        for (int iEm = 1; iEm <= numEm_; ++iEm) sfval_[idx(iEm, iPm)] = sfval_[idx(iEm, iPm)] / sftotnorm_;

    // Corrected mode (FORTRAN-BUG D-033 F3): the differential spectral function at every node, with the node's own
    // bin widths, S = SF_PROB / (4 pi P^2 dP dE) (the definition of SF_PROB in the sf_lookup.f header, l.8-11).
    dens_.assign(sfval_.size(), 0.0);
    for (int iPm = 1; iPm <= numPm_; ++iPm)
        for (int iEm = 1; iEm <= numEm_; ++iEm)
            dens_[idx(iEm, iPm)] =
                sfval_[idx(iEm, iPm)] / (4. * constants::pi * Pmval_[iPm] * Pmval_[iPm] * dPm_[iPm] * dEm_[iEm]);
    Em1_ = Em2_ = sf1_ = sf2_ = 0.0;
    gx_.assign(nEmmax + 1, 0.0);
    gy_.assign(nEmmax + 1, 0.0);
    return true;
}

// =====================================================================================================================
// Fortran: sf_lookup.f:sf_lookup (l.99-190), applied to table t
// =====================================================================================================================
double SfTable::interpolate(const std::vector<double>& t, double Em, double Pm, bool stale) {
    if (numPm_ < 2) throw std::logic_error("sf_lookup: the spectral-function table is not loaded (sf_lookup_init)");
    int iPm;
    double w1, w2;
    // find nearest Pm value (l.118-154)
    if (Pm >= Pmval_[numPm_]) {
        iPm = numPm_ - 1;
        w1 = 0;
        w2 = 1;
    } else if (Pm <= Pmval_[1]) {
        iPm = 1;
        w1 = 1;
        w2 = 0;
    } else {
        int ind = 1;
        while (Pm > Pmval_[ind]) ind = ind + 1;
        iPm = ind - 1;
        w2 = (Pm - Pmval_[iPm]) / (Pmval_[iPm + 1] - Pmval_[iPm]);
        w1 = (Pmval_[iPm + 1] - Pm) / (Pmval_[iPm + 1] - Pmval_[iPm]);
        if (std::abs(w1 * Pmval_[iPm] + w2 * Pmval_[iPm + 1] - Pm) > 0.0001)
            throw FortranStop("sf_lookup: w1,w2,Pm,Pmval(iPm) inconsistent (sf_lookup.f:135-138)");
    }
    if (std::abs(w1 + w2 - 1) > 0.0001) throw FortranStop("sf_lookup: iPm,Pm,w1+w2 (sf_lookup.f:155-158)");

    // l.160-182. SAVEd Em1, Em2, sf1, sf2 (fortran_compat) or fresh locals (corrected).
    double Em1 = stale ? Em1_ : 0., Em2 = stale ? Em2_ : 0., sf1 = stale ? sf1_ : 0., sf2 = stale ? sf2_ : 0.;
    bool assigned = false;
    if (Em <= Emval_[1]) {   // linear extrapolation
        Em1 = Emval_[1];
        Em2 = Emval_[2];
        sf1 = w1 * t[idx(1, iPm)] + w2 * t[idx(1, iPm + 1)];
        sf2 = w1 * t[idx(2, iPm)] + w2 * t[idx(2, iPm + 1)];
        assigned = true;
    } else if (Em > Emval_[numEm_]) {   // linear extrapolation
        Em1 = Emval_[numEm_ - 1];
        Em2 = Emval_[numEm_];
        sf1 = w1 * t[idx(numEm_ - 1, iPm)] + w2 * t[idx(numEm_ - 1, iPm + 1)];
        sf2 = w1 * t[idx(numEm_, iPm)] + w2 * t[idx(numEm_, iPm + 1)];
        assigned = true;
    } else {
        for (int iEm = 1; iEm <= numEm_ - 1; ++iEm) {
            // FORTRAN-BUG(D-033 F6): original Fortran sf_lookup.f:174: `if (Em.ge.Emval(iEm) .and. Em.lt.Emval(iEm+1))`
            // -- Em == Emval(numEm) matches no interval and the SAVEd Em1/Em2/sf1/sf2 of the previous call are used
            // (NaN on the first call); corrected because the last node belongs to the last interval.
            const bool inLast = !stale && iEm + 1 == numEm_ && Em == Emval_[numEm_];
            if ((Em >= Emval_[iEm] && Em < Emval_[iEm + 1]) || inLast) {
                Em1 = Emval_[iEm];
                Em2 = Emval_[iEm + 1];
                sf1 = w1 * t[idx(iEm, iPm)] + w2 * t[idx(iEm, iPm + 1)];
                sf2 = w1 * t[idx(iEm + 1, iPm)] + w2 * t[idx(iEm + 1, iPm + 1)];
                assigned = true;
            }
        }
    }
    if (stale) {
        Em1_ = Em1;
        Em2_ = Em2;
        sf1_ = sf1;
        sf2_ = sf2;
    } else if (!assigned) {
        return 0.0;   // corrected: NaN Em has no table value (no stale state)
    }
    // FORTRAN-QUIRK(D-033 F6): called logsf, but the interpolation is linear (no log is taken).
    const double logsf = (sf1 + (Em - Em1) * (sf2 - sf1) / (Em2 - Em1));
    double SF = logsf;
    if (SF < 1.e-20) SF = 0;
    return SF;
}

double SfTable::lookup(double Em, double Pm) {
    if (!compat_ && (std::isnan(Em) || std::isnan(Pm))) return 0.0;
    return interpolate(sfval_, Em, Pm, compat_);
}

// Fortran: sf_lookup.f:sf_lookup_diff (l.85-95)
double SfTable::lookupDiff(double Em, double Pm) {
    if (compat_) {
        // FORTRAN-BUG(D-033 F3): original Fortran sf_lookup.f:93: `SFd = SF/4/3.1415926535/(Pm*Pm)/5.0/20.0`; corrected
        // (default mode) because the bin probability must be divided by the volume of ITS bin (4 pi Pmval^2 dPm dEm
        // with the widths of the file), not by the event's own Pm^2 (+Inf/NaN at Pm = 0) and fixed 5 x 20 widths.
        const double SF = lookup(Em, Pm);
        return SF / 4 / 3.1415926535 / (Pm * Pm) / 5.0 / 20.0;
    }
    if (std::isnan(Em) || std::isnan(Pm)) return 0.0;
    return interpolate(dens_, Em, Pm, false);
}

// =====================================================================================================================
// Fortran: sf_lookup.f:generate_em (l.195-245)
// =====================================================================================================================
double SfTable::generateEm(double Pm, IRandom& rng) {
    std::vector<double>& x = gx_;
    std::vector<double>& y = gy_;
    x[1] = Emval_[1];
    y[1] = lookup(x[1], Pm);
    for (int iEm = 2; iEm <= numEm_; ++iEm) {
        x[iEm] = Emval_[iEm];
        y[iEm] = lookup(x[iEm], Pm);
        y[iEm] = y[iEm] + y[iEm - 1];
    }
    // Normalize the Em distribution. (y(numEm) itself becomes 1 at the last iteration: later entries use it.)
    for (int iEm = 1; iEm <= numEm_; ++iEm) y[iEm] = y[iEm] / y[numEm_];
    const double ranprob = rng.uniform();
    int ind = 1;
    while (ranprob > y[ind]) {
        ind = ind + 1;
        if (ind > nEmmax) throw std::logic_error("generate_em: Em search ran past the table");
    }
    return Emval_[ind] + dEm_[ind] * (rng.uniform() - 0.5);
}

// =====================================================================================================================
// Fortran: event.f:complete_main (l.1402-1430)
// =====================================================================================================================
double TheorySfWeight(const TheoryData& t, bool doing_heavy, double Em, double Pm) {
    using constants::pi;
    double SF_weight = 0.0;
    for (int i = 0; i < t.nrhoPm; ++i) {   // Fortran i = 1..nrhoPm -> i+1
        double weight = 0.0;
        const TheoryData::PmAxis& ax = t.Pm[static_cast<std::size_t>(i)];
        // ... use linear interpolation to determine rho(pm)
        const double r = (Pm - ax.min) / ax.bin;
        if (r >= 0 && r <= ax.n) {
            int iPm1 = static_cast<int>(std::lround(r));   // gfortran nint: half away from zero
            if (iPm1 == 0) iPm1 = 1;
            if (iPm1 == ax.n) iPm1 = ax.n - 1;
            const double frac = r + 0.5 - static_cast<double>(iPm1);
            const double b = t.theory[static_cast<std::size_t>(i)][static_cast<std::size_t>(iPm1)];
            const double a = t.theory[static_cast<std::size_t>(i)][static_cast<std::size_t>(iPm1 + 1)] - b;
            weight = a * frac + b;
        }
        if (doing_heavy) {
            const double width = t.Emsig[static_cast<std::size_t>(i)] / 2.0;
            if (Em < t.E_Fermi) weight = 0.0;
            const double d = Em - t.Em[static_cast<std::size_t>(i)];
            weight = weight / pi / t.Em_int[static_cast<std::size_t>(i)] * width / (d * d + width * width);
        }
        SF_weight = SF_weight + weight * t.nprot[static_cast<std::size_t>(i)];
    }
    return SF_weight;
}

}  // namespace nuclear
}  // namespace simc
