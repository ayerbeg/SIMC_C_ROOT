// src/setup/SetupTheory.cpp — part of BuildRunSetup (see SetupInternal.h): the SIMC data-file readers
// (init.f:theory_init, sf_lookup_init's Emval(numEm), the momentum distributions) and BuildRunSetup phase 2:
// init.f:target_init through the target model's ITargetLossProvider, theory_init, the momentum distribution and the
// spectral-function choice (dbase.f:551-681).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <sstream>

#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/core/Units.h"
#include "SetupInternal.h"

namespace simc {

namespace {

namespace K = constants;
using setup_detail::nint;

// ---------------------------------------------------------------------------------------------
// List-directed record reader (just enough of Fortran `read(u,*)` for the SIMC data files):
// every read starts on a new record, takes values across records if needed, and discards the
// rest of the last record it touched.
// ---------------------------------------------------------------------------------------------
class ListReader {
public:
    explicit ListReader(std::istream& in) {
        std::string l;
        while (std::getline(in, l)) lines_.push_back(l);
    }
    enum Status { kOk, kEnd, kErr };
    Status read(std::vector<double>& v, std::size_t n) {
        v.clear();
        while (v.size() < n) {
            if (next_ >= lines_.size()) return kEnd;
            std::string l = lines_[next_++];
            for (char& c : l)
                if (c == ',') c = ' ';
            std::istringstream is(l);
            std::string tok;
            while (v.size() < n && is >> tok) {
                // Fortran D exponents
                for (char& c : tok)
                    if (c == 'd' || c == 'D') c = 'E';
                char* end = nullptr;
                const double x = std::strtod(tok.c_str(), &end);
                if (end == tok.c_str() || *end != '\0') return kErr;
                v.push_back(x);
            }
        }
        return kOk;
    }

private:
    std::vector<std::string> lines_;
    std::size_t next_{0};
};

}  // namespace


// =============================================================================================
// Fortran: init.f:theory_init (l.828-920)
// =============================================================================================
bool ReadTheoryFile(const std::string& path, bool doing_heavy, TheoryData& t, std::vector<std::string>* msgs) {
    constexpr int ntheorymax = 500;   // simulate.inc:117
    t = TheoryData{};
    t.file = path;
    std::ifstream f(path);
    if (!f) {
        if (msgs) msgs->push_back("ERROR: theory_init cannot open " + path);
        return false;
    }
    ListReader rd(f);
    std::vector<double> v;
    auto fail = [&]() {
        if (msgs) msgs->push_back("ERROR: theory_init failed to read in the theory file");
        return false;
    };
    // read(1,*,err=40) nrhoPm, absorption, E_Fermi
    if (rd.read(v, 3) != ListReader::kOk) return fail();
    t.nrhoPm = static_cast<int>(v[0]);
    const double absorption = v[1];
    t.E_Fermi = v[2];
    if (t.nrhoPm < 0) return fail();
    t.nprot.resize(t.nrhoPm);
    t.Em.resize(t.nrhoPm);
    t.Emsig.resize(t.nrhoPm);
    t.bs_norm.resize(t.nrhoPm);
    t.Pm.resize(t.nrhoPm);
    t.theory.assign(t.nrhoPm + 1, std::vector<double>(ntheorymax + 1, 0.0));   // 1-based n
    for (int m = 0; m < t.nrhoPm; ++m) {
        if (rd.read(v, 4) != ListReader::kOk) return fail();
        t.nprot[m] = v[0] * absorption;
        t.Em[m] = v[1];
        t.Emsig[m] = v[2];
        t.bs_norm[m] = v[3];
    }
    std::vector<double> Pm_values(ntheorymax + 2, 0.0);   // 1-based, static (zero) in the Fortran
    if (rd.read(v, 2) != ListReader::kOk) return fail();   // err=40 (no end= : EOF would abort)
    Pm_values[1] = v[0];
    t.theory[0][1] = v[1];
    for (int m = 0; m < t.nrhoPm; ++m) {
        int n = 2;
        auto st = rd.read(v, 2);
        if (st == ListReader::kErr) return fail();
        if (st == ListReader::kOk) {
            Pm_values[2] = v[0];
            t.theory[m][2] = v[1];
            while (Pm_values[n] > Pm_values[n - 1]) {
                n = n + 1;
                if (n > ntheorymax) return fail();   // Fortran would overrun theory(,ntheorymax)
                st = rd.read(v, 2);
                if (st == ListReader::kErr) return fail();
                if (st == ListReader::kEnd) break;   // end=50
                Pm_values[n] = v[0];
                t.theory[m][n] = v[1];
            }
        }
        // 50: figure out details of the Pm axes
        TheoryData::PmAxis& ax = t.Pm[m];
        ax.n = n - 1;
        ax.bin = Pm_values[2] - Pm_values[1];
        ax.min = Pm_values[1] - ax.bin / 2.;
        ax.max = Pm_values[ax.n] + ax.bin / 2.;
        if (std::abs(ax.min + ax.bin * ax.n - ax.max) > 0.1) {
            if (msgs) msgs->push_back("ERROR: theory_init found unequal Pm bins in distribution number " +
                                      std::to_string(m + 1) + "!");
            return false;
        }
        // prepare for the next loop
        Pm_values[1] = Pm_values[n];
        t.theory[m + 1][1] = t.theory[m][n];
    }
    // renormalize the momentum distributions
    for (int m = 0; m < t.nrhoPm; ++m)
        for (int n = 1; n <= t.Pm[m].n; ++n) t.theory[m][n] = t.theory[m][n] / t.bs_norm[m];
    // integral of the Em distribution above E_Fermi
    t.Em_int.assign(t.nrhoPm, 0.0);
    if (doing_heavy) {
        for (int m = 0; m < t.nrhoPm; ++m)
            t.Em_int[m] = (K::pi / 2. + std::atan((t.Em[m] - t.E_Fermi) / (0.5 * t.Emsig[m]))) / K::pi;
    }
    t.loaded = true;
    return true;
}

// =============================================================================================
// Fortran: sf_lookup.f:sf_lookup_init (l.22-49), the part the set-up needs: Emval(numEm), the Em of the last line of
// the first Pm block (iPm = 1). The table itself is read by nuclear::SfTable (physics library) for the reaction.
// =============================================================================================
bool ReadBenharEmMax(const std::string& path, double& emmax) {
    std::ifstream f(path);
    if (!f) return false;
    ListReader rd(f);
    std::vector<double> v;
    if (rd.read(v, 2) != ListReader::kOk) return false;
    const int numEm = static_cast<int>(v[1]);
    if (numEm < 1) return false;
    for (int iEm = 1; iEm <= numEm; ++iEm) {
        if (rd.read(v, 6) != ListReader::kOk) return false;
        emmax = v[1];   // Emval(iEm) = tmpEm (iPm = 1)
    }
    return true;
}

// =============================================================================================
// Fortran: dbase.f:578-586 (momentum distribution pval/mprob)
// =============================================================================================
bool ReadMomentumDistribution(const std::string& path, MomentumDistribution& d) {
    d = MomentumDistribution{};
    d.file = path;
    std::ifstream f(path);
    if (!f) return false;   // Fortran: open(status='old') aborts
    ListReader rd(f);
    std::vector<double> v;
    for (int ii = 1; ii <= 2000; ++ii) {
        const auto st = rd.read(v, 2);
        if (st != ListReader::kOk) break;   // end=999 (an error would abort the Fortran)
        d.pval.push_back(v[0]);
        d.mprob.push_back(v[1]);
    }
    if (d.mprob.empty()) return false;
    const double last = d.mprob.back();
    for (double& x : d.mprob) x = x / last;
    // NB: mprob(nump) itself becomes exactly 1 on the last iteration; earlier entries used the
    // original last value, as in the Fortran loop (mprob(nump) is only overwritten at ii=nump).
    return true;
}

namespace setup_detail {

bool BuildTargetAndTheory(RunSetup& s, const SimcParameters& p, const SetupOptions& opt, ITargetLossProvider* losses,
                          double& benharEmMax) {
    auto stop = [&](const std::string& why) { s.stops.push_back("STOP " + why); };
    auto msg = [&](const std::string& m) { s.messages.push_back(m); };
    const ReactionFlags& f = s.flags;
    RunSwitches& sw = s.sw;
    TargetInfo& t = s.targ;
    const int A = nint(p.targ.A);
    {
        TargetInitArgs a;
        a.electron_arm = p.electron_arm;
        a.hadron_arm = p.hadron_arm;
        a.using_Eloss = sw.using_Eloss;
        a.using_Coulomb = sw.using_Coulomb;
        a.the_cent = s.spec.e.theta;
        a.thp_cent = s.spec.p.theta;
        a.Pp_cent = s.spec.p.P;
        a.Mh2 = s.Mh2;
        a.Ebeam = s.Ebeam;
        a.Pe_cent = s.spec.e.P;
        losses->targetInit(a, t, s.messages);
    }

    s.data_dir = opt.data_dir;
    // ---- dbase.f:555-558 theory_init
    if (f.doing_deuterium || (f.doing_heavy && !p.use_benhar_sf)) {
        // init.f:839-851. D-033 F1: every A other than 2/12/56/197 silently gets c12.theory in the Fortran. SIMC++ keeps
        // that default but says so, and lets the input name the file (SetupOptions::sf_file).
        std::string file;
        SpectralFunctionChoice& sf = s.sf;
        sf.type = "theory";
        if (!opt.sf_file.empty()) {
            file = opt.sf_file;
            sf.explicit_file = true;
        } else if (A == 2) file = "h2.theory";
        else if (A == 12) file = "c12.theory";
        else if (A == 56) file = "fe56.theory";
        else if (A == 197) file = "au197.theory";
        else {
            msg("No Theory File (spectral function) for A = " + std::to_string(t.A) + ", defaulting to c12.theory");
            msg("WARNING: the spectral function of this target is the C-12 one (c12.theory: 6 protons in the C-12 shells, "
                "E_Fermi 15.96 MeV), as in the Fortran (init.f:839-851). Name the file explicitly with "
                "pipeline.reaction.spectral_function.file to make this choice deliberate.");
            file = "c12.theory";
            sf.fallback = true;
        }
        sf.file = file;
        const bool okT = ReadTheoryFile(joinPath(opt.data_dir, file), f.doing_heavy, s.theory, &s.messages);
        s.theory.file = file;
        if (!okT) {
            stop("THEORY_INIT failed!");
            return false;
        }
        // init.f:891 FORTRAN: theory_init overwrites doing_deuterium from the file contents
        s.flags.doing_deuterium = s.theory.nrhoPm == 1 && s.theory.E_Fermi < 1.0;
    }

    // ---- dbase.f:562-587 momentum distribution
    // FORTRAN-BUG(D-048 F1): original Fortran dbase.f:563-565 and 594 (momentum distribution and spectral function):
    // `if(doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon.or.doing_deutsemi)` /
    // `if(doing_hepi.or.doing_hekaon .or. ...)` omit doing_deutdelta / doing_hedelta, although generate (event.f:336-371)
    // draws the D/A(e,e'p)pi Fermi momentum from pval/mprob and the A > 2 missing energy with generate_em (sf_lookup),
    // and limits_init uses pval(nump) / Emval(numEm) (init.f:349-358). The Fortran never reaches them (the D/A Delta runs
    // STOP, dbase.f:779-788); corrected together with that STOP: the delta flags load the same files as the pion.
    const bool deutdelta = f.doing_deutdelta && !s.fortran_compat, hedelta = f.doing_hedelta && !s.fortran_compat;
    if (f.doing_deutpi || f.doing_hepi || f.doing_deutkaon || f.doing_hekaon || f.doing_deutsemi || deutdelta || hedelta) {
        std::string file;
        if (f.doing_deutpi || f.doing_deutkaon || f.doing_deutsemi || deutdelta) file = "deut.dat";
        if (f.doing_hepi || f.doing_hekaon || hedelta) {
            if (A == 3) file = "he3.dat";
            else if (A == 4) file = "he4.dat";
            else if (A == 12) file = "c12.dat";
            else {
                msg("No Momentum Distribution for A = " + std::to_string(t.A) + ", defaulting to carbon");
                file = "c12.dat";
            }
        }
        if (!ReadMomentumDistribution(joinPath(opt.data_dir, file), s.pdist)) {
            stop("cannot read momentum distribution " + file);
            return false;
        }
        s.pdist.file = file;
    }

    // ---- dbase.f:592-621 spectral function (sf_lookup_init). For A(e,e'p) with use_benhar_sf the choice of the file and
    // of the proton/neutron column is made here; the table is read by the "benhar" ISpectralFunction (sf_lookup_init).
    // (also for A > 2 quasifree pion/kaon production, doing_hepi/doing_hekaon: the table of generate_em)
    benharEmMax = 0.0;   // Emval(numEm)
    if ((f.doing_heavy && p.use_benhar_sf) || f.doing_hepi || f.doing_hekaon || hedelta) {
        SpectralFunctionChoice& sf = s.sf;
        sf.type = "benhar";
        if (!opt.sf_file.empty()) {
            sf.file = opt.sf_file;
            sf.explicit_file = true;
        } else if (A == 3) {
            msg("Using the mod version of 3He S.F. rather than Paris.");
            sf.file = "benharsf_3mod.dat";
        } else if (A == 4) sf.file = "benharsf_4.dat";
        else if (A == 12) sf.file = "benharsf_12.dat";
        else if (A == 56) sf.file = "benharsf_56.dat";
        else if (A == 197) sf.file = "benharsf_197.dat";
        else {
            msg("No spectral function for A = " + std::to_string(t.A) + ", Defaulting to Carbon S.F.");
            msg("WARNING: the spectral function of this target is the C-12 Benhar one (benharsf_12.dat), as in the "
                "Fortran (dbase.f:611-614). Name the file explicitly with pipeline.reaction.spectral_function.file.");
            sf.file = "benharsf_12.dat";
            sf.fallback = true;
        }
        if (std::abs(t.Mtar_struck - K::Mp) <= 1.e-6) {
            sf.proton = true;
        } else if (std::abs(t.Mtar_struck - K::Mn) <= 1.e-6) {
            sf.proton = false;
        } else {
            msg("targ%Mtar_struck not equal to Mp or Mn: DEFAULTING TO PROTON SPECTRAL FUNCTION!!!");
            sf.proton = true;
        }
    }
    if (f.doing_hepi || f.doing_hekaon || hedelta) {   // sf_lookup_init: Emval(numEm) for VERTEXedge%Em%max (init.f:356)
        if (!ReadBenharEmMax(joinPath(opt.data_dir, s.sf.file), benharEmMax)) {
            stop("cannot read spectral function " + s.sf.file);
            return false;
        }
    }
    // ---- dbase.f:624-681 kaon model tables: they belong to the kaon cross-section model
    // (doing_kaon: read by the "kaon" physics model, kaon::SaghaiTables)
    return true;
}

}  // namespace setup_detail
}  // namespace simc
