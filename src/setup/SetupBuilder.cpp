// src/setup/SetupBuilder.cpp
//
// Port of the Fortran start-up sequence (simc.f l.53-101 without calculate_central):
// dbase.f:dbase_read post-read logic, init.f:target_init (non energy-loss part), theory_init,
// limits_init, radc_init, simc.f:94-101. Line numbers refer to simc_gfortran/{dbase,init,simc}.f.
//
// Numerics: the Fortran is compiled with -fdefault-real-8, so every literal is a double and
// every expression below keeps the Fortran operation order (left to right as written).
#include "simc/setup/SetupBuilder.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "simc/core/SimcConstants.h"

namespace simc {

namespace {

namespace K = constants;

int nint(double x) { return static_cast<int>(std::lround(x)); }

std::string joinPath(const std::string& dir, const std::string& file) {
    if (dir.empty()) return file;
    return dir.back() == '/' ? dir + file : dir + "/" + file;
}

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

// =============================================================================================
// BuildRunSetup
// =============================================================================================
RunSetup BuildRunSetup(const SimcParameters& p, const SetupOptions& opt, ITargetLossProvider* losses) {
    RunSetup s = RunSetup::FromParameters(p, opt.run_name);
    auto stop = [&](const std::string& why) { s.stops.push_back("STOP " + why); };
    auto msg = [&](const std::string& m) { s.messages.push_back(m); };
    const int A = nint(p.targ.A);
    double mrec_guess = 0.0;   // dbase.f local 'Mrec_guess' (static zero under -fno-automatic)

    // ---- dbase.f:124-222 flags (DeriveReactionFlags) and Mh. Mh is chosen BEFORE doing_phsp clears
    // the reaction flags (dbase.f:227-234), so evaluate the branch on the pre-phsp flags.
    std::vector<std::string> flag_stops;
    s.flags = DeriveReactionFlags(p, &flag_stops);
    if (!flag_stops.empty()) {
        s.stops = flag_stops;
        return s;
    }
    SimcParameters pre = p;
    pre.doing_phsp = false;
    const ReactionFlags f0 = DeriveReactionFlags(pre);
    if (f0.doing_pion) {                         // dbase.f:138-139
        s.Mh = K::Mpi;
        if (f0.doing_pizero) s.Mh = K::Mpi0;
    } else if (f0.doing_kaon) {                  // dbase.f:163
        s.Mh = K::Mk;
    } else if (f0.doing_delta) {                 // dbase.f:182-185
        s.Mh = K::Mp;
        if (A >= 2) msg("WARNING: Delta cross section model only set up for proton target!");
    } else if (f0.doing_semi) {                  // dbase.f:191-205
        if (f0.doing_semipi) {
            s.Mh = K::Mpi;
            if (f0.doing_pizero) s.Mh = K::Mpi0;
        } else if (f0.doing_semika) {
            s.Mh = K::Mk;
        }
        if (f0.doing_hydsemi && p.do_fermi) msg("WARNING: Cannot do Fermi motion for Hydrogen!");
    } else if (f0.doing_rho) {                   // dbase.f:208
        s.Mh = K::Mrho;
    } else {                                     // dbase.f:216
        s.Mh = K::Mp;
    }
    s.Mh2 = s.Mh * s.Mh;                         // dbase.f:224-225
    s.Mh2_final = s.Mh2;
    const ReactionFlags& f = s.flags;            // post-phsp flags from here on

    // ---- switches that later code needs (convert_to_logical values)
    RunSwitches& sw = s.sw;
    sw.using_Eloss = p.using_Eloss;
    sw.correct_Eloss = p.correct_Eloss;
    sw.using_Coulomb = p.using_Coulomb;
    sw.mc_smear = p.mc_smear;
    sw.using_rad = p.using_rad;
    sw.rad_flag = p.rad_flag;
    sw.extrad_flag = p.extrad_flag;
    if (p.doing_phsp) sw.rad_flag = 0;           // dbase.f:227-228

    // ---- dbase.f:238-263 kinematics_main
    s.Ebeam = p.Ebeam;
    s.dEbeam = p.Ebeam * p.dEbeam / 100.;
    auto arm = [&](const SimcParameters::Arm& a, SpecArm& o) {
        o.P = a.P;
        o.theta = std::abs(a.theta) / K::degrad;
        o.cos_th = std::cos(o.theta);
        o.sin_th = std::sin(o.theta);
        o.offset.x = a.offset.x;
        o.offset.y = a.offset.y;
        o.offset.z = a.offset.z;
        o.offset.xptar = a.offset.xptar;   // converted below (dbase.f:472-475)
        o.offset.yptar = a.offset.yptar;
    };
    arm(p.spec.e, s.spec.e);
    arm(p.spec.p, s.spec.p);
    auto armPhi = [&](int code, double& phi, const char* which) {
        if (code == 1 || code == 3 || code == 7) {
            phi = 3 * K::pi / 2.;
        } else if (code == 2 || code == 4 || code == 5 || code == 6 || code == 8) {
            phi = K::pi / 2.;
        } else {
            stop(std::string("I dont know what phi should be for the ") + which + " arm");
            return false;
        }
        return true;
    };
    if (!armPhi(p.electron_arm, s.spec.e.phi, "electron")) return s;   // dbase.f:247-254
    if (!armPhi(p.hadron_arm, s.spec.p.phi, "hadron")) return s;       // dbase.f:256-263

    // ---- dbase.f:267-285 target masses
    TargetInfo& t = s.targ;
    t.A = p.targ.A;
    t.Z = p.targ.Z;
    t.mass_amu = p.targ.mass_amu;
    t.mrec_amu = p.targ.mrec_amu;
    t.rho = p.targ.rho;
    t.thick = p.targ.thick;
    t.angle = p.targ.angle;
    t.abundancy = p.targ.abundancy;
    t.zoffset = p.targ.zoffset;
    t.xoffset = p.targ.xoffset;
    t.yoffset = p.targ.yoffset;
    t.fr1 = p.targ.fr1;
    t.fr2 = p.targ.fr2;
    t.fr_pattern = p.targ.fr_pattern;
    t.can = p.targ.can;
    t.N = t.A - t.Z;
    t.M = t.mass_amu * K::amu;
    t.Mrec = t.mrec_amu * K::amu;
    if (A == 1) {
        if (std::abs(t.M - K::Mp) > 0.1) msg("WARNING: forcing target mass to Mp!!!");
        if (std::abs(t.Mrec) > 0.1) msg("WARNING: forcing recoil mass to zero!!!");
        t.M = K::Mp;
        t.Mrec = 0.;
    } else if (A == 2) {
        if (std::abs(t.M - K::Md) > 0.1) msg("WARNING: forcing target mass to Md!!!");
        t.M = K::Md;
    } else {
        mrec_guess = t.M - K::Mp;   // approx A-1 mass
        if (std::abs(t.Mrec - mrec_guess) > 100.) {
            msg("targ.Mrec replaced by the expectation M - Mp (dbase.f:283-284)");
            t.Mrec = mrec_guess;
        }
    }

    // ---- dbase.f:290-439 struck / recoil masses and sign_hadron (post-phsp flags:
    // FORTRAN-QUIRK(D-009): with doing_phsp no branch is taken, Mtar_struck = Mrec_struck = 0).
    if (f.doing_eep) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = 0.0;
        s.sign_hadron = 1.0;
    } else if (f.doing_delta) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mpi;
        s.sign_hadron = 1.0;
    } else if (f.doing_semi) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mp;
        s.sign_hadron = p.doing_hplus ? 1.0 : -1.0;
    } else if (f.doing_rho) {
        t.Mtar_struck = K::Mp;
        t.Mrec_struck = K::Mp;
        s.sign_hadron = p.doing_hplus ? 1.0 : -1.0;
    } else if (f.doing_pion) {
        const int w = p.which_pion;
        if (f.doing_pizero) {
            if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mp; s.sign_hadron = 1.0; }
            else if (w == 1) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Mn; s.sign_hadron = -1.0; }
            else if (w == 2) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = 1.0; }
            else if (w == 3) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::MDelta; s.sign_hadron = -1.0; }
        } else {
            if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mn; s.sign_hadron = 1.0; }
            else if (w == 1) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Mp; s.sign_hadron = -1.0; }
            else if (w == 2) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = 1.0; }
            else if (w == 3) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::MDelta; s.sign_hadron = -1.0; }
            else if (w == 10) {
                t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec;
                mrec_guess = t.M - K::Mp + K::Mn; s.sign_hadron = 1.0;
            } else if (w == 11) {
                t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec;
                mrec_guess = t.M - K::Mn + K::Mp; s.sign_hadron = -1.0;
            } else {
                stop("Bad value for which_pion");
                return s;
            }
        }
        if (w >= 10) {   // coherent pion production (dbase.f:382-391)
            if (std::abs(t.Mrec_struck - mrec_guess) > 100.) {
                t.Mrec_struck = mrec_guess;
                msg("For coherent pion production the recoil mass was replaced by M_A + M_N(struck) - M_N(recoil)");
            }
            t.Mrec = 0.;
        }
    } else if (f.doing_kaon) {
        const int w = p.which_kaon;
        if (w == 0) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Mlambda; s.sign_hadron = 1.0; }
        else if (w == 1) { t.Mtar_struck = K::Mp; t.Mrec_struck = K::Msigma0; s.sign_hadron = 1.0; }
        else if (w == 2) { t.Mtar_struck = K::Mn; t.Mrec_struck = K::Msigma_minus; s.sign_hadron = 1.0; }
        else if (w == 10) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mp + K::Mlambda; s.sign_hadron = 1.0;
        } else if (w == 11) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mp + K::Msigma0; s.sign_hadron = 1.0;
        } else if (w == 12) {
            t.Mtar_struck = t.M; t.Mrec_struck = t.Mrec; mrec_guess = t.M - K::Mn + K::Msigma_minus; s.sign_hadron = 1.0;
        } else {
            stop("Bad value for which_kaon");
            return s;
        }
        if (w >= 10) {   // dbase.f:427-436
            if (std::abs(t.Mrec_struck - mrec_guess) > 100.) {
                t.Mrec_struck = mrec_guess;
                msg("For coherent kaon production the recoil mass was replaced by M_A - M_N + M_Y");
            }
            t.Mrec = 0.;
        }
    }

    // ---- dbase.f:443-448 deuterium recoil mass
    if (A == 2) {
        mrec_guess = K::Mp + K::Mn - t.Mtar_struck;
        if (std::abs(t.Mrec - mrec_guess) > 0.1) msg("WARNING: forcing recoil mass to M_Nucleon");
        t.Mrec = mrec_guess;
    }

    // ---- dbase.f:451-468 target geometry
    t.thick = t.thick / 1000.;
    t.length = t.thick / t.rho;
    t.angle = t.angle / K::degrad;
    s.targ_Bangle = p.targ_Bangle / K::degrad;
    s.targ_Bphi = p.targ_Bphi / K::degrad;
    if (t.Z < 2.4) {   // liquid target
        if (std::abs(t.angle) > 0.0001) {
            t.angle = 0.0;
            msg("Forcing target angle to zero for cryotarget.");
        }
        if (t.can != 1 && t.can != 2 && t.can != 3) {
            stop("bad targ.can value");
            return s;
        }
    }
    if (std::sin(t.angle) > 0.85) {
        stop("BAD targ.angle (0 is perp. to beam, +ve is rotated towards SOS)");
        return s;
    }

    // ---- dbase.f:472-475 offsets mrad -> rad
    s.spec.e.offset.xptar = s.spec.e.offset.xptar / 1000.;
    s.spec.e.offset.yptar = s.spec.e.offset.yptar / 1000.;
    s.spec.p.offset.xptar = s.spec.p.offset.xptar / 1000.;
    s.spec.p.offset.yptar = s.spec.p.offset.yptar / 1000.;

    // ---- dbase.f:479-488 SPedge (min_max_init ±1e10 defaults already in params)
    auto spedge = [](const SimcParameters::ArmEdge& in, ArmCuts& o) {
        o.delta.min = in.delta.min;
        o.delta.max = in.delta.max;
        if (o.delta.min <= -100.0) o.delta.min = -99.99;
        o.yptar.min = in.yptar.min / 1000.;
        o.yptar.max = in.yptar.max / 1000.;
        o.xptar.min = in.xptar.min / 1000.;
        o.xptar.max = in.xptar.max / 1000.;
    };
    spedge(p.SPedge.e, s.SPedge.e);
    spedge(p.SPedge.p, s.SPedge.p);

    // ---- dbase.f:494-518 radiation tails
    const int ot = p.one_tail;
    sw.doing_tail[0] = ot == 0 || ot == 1 || ot == -2 || ot == -3;
    sw.doing_tail[1] = ot == 0 || ot == 2 || ot == -3 || ot == -1;
    sw.doing_tail[2] = ot == 0 || ot == 3 || ot == -1 || ot == -2;
    if (!sw.doing_tail[0] || !sw.doing_tail[1]) msg("**WARNING: Shutting off tails 1 or 2 not well defined");
    if (std::abs(ot) > 3 && p.using_rad) {
        stop("Moron! one_tail>3 turns radiation off, but using_rad wants it on.");
        return s;
    }
    if (!p.using_rad) sw.doing_tail = {false, false, false};
    sw.hardwired_rad = p.Egamma_gen_max > 0.01;

    // ---- dbase.f:521-528 gen angles /1000 (only the ±1e10 min_max_init values at this point)
    GenLimits& gen = s.gen;
    gen.xwid = p.gen.xwid;
    gen.ywid = p.gen.ywid;
    gen.e.yptar.min = gen.e.yptar.min / 1000.;
    gen.e.yptar.max = gen.e.yptar.max / 1000.;
    gen.p.yptar.min = gen.p.yptar.min / 1000.;
    gen.p.yptar.max = gen.p.yptar.max / 1000.;
    gen.e.xptar.min = gen.e.xptar.min / 1000.;
    gen.e.xptar.max = gen.e.xptar.max / 1000.;
    gen.p.xptar.min = gen.p.xptar.min / 1000.;
    gen.p.xptar.max = gen.p.xptar.max / 1000.;

    // ---- dbase.f:531-547 flags
    sw.using_E_arm_montecarlo = p.spect_mode != 1 && p.spect_mode != -1;
    sw.using_P_arm_montecarlo = p.spect_mode != 1 && p.spect_mode != -2;
    s.cuts.Em.min = p.cuts.Em.min;
    s.cuts.Em.max = p.cuts.Em.max;
    if (f.doing_pion || f.doing_kaon || f.doing_delta || (s.cuts.Em.min == s.cuts.Em.max)) {
        s.cuts.Em.min = -1.e6;
        s.cuts.Em.max = 1.e6;
    }
    if (std::abs(p.deForest_flag) > 1) {
        stop("Idiot! check setting of deForest_flag");
        return s;
    }
    if (sw.correct_Eloss && !sw.using_Eloss) {
        msg("using_Eloss is false, SO WE ARE FORCING CORRECT_ELOSS TO BE FALSE!!!");
        sw.correct_Eloss = false;
    }
    if (nint(t.Z) == 1) sw.using_Coulomb = false;   // no coulomb corr. for Z=1

    // ---- init.f:target_init (l.1-88), called at dbase.f:551
    t.L1 = std::log(184.15) - std::log(t.Z) / 3.0;
    t.L2 = std::log(1194.) - 2. * std::log(t.Z) / 3.0;
    if (t.Z == 1) {
        t.L1 = 5.31;
        t.L2 = 6.144;
    }
    {
        const double za2 = (t.Z * K::alpha) * (t.Z * K::alpha);
        const double fc = za2 * (1.202 + za2 * (-1.0369 + za2 * 1.008 / (za2 + 1)));
        if (A == 1) {
            t.X0 = 61.28;
        } else if (A == 2) {
            t.X0 = 122.4;
        } else if (A == 4) {
            t.X0 = 94.32;
        } else {
            if (A == 3) msg("Using Tsai formula for He-3 radiation length!!! (is there a better value?)");
            t.X0 = 716.405 * t.A / t.Z / (t.Z * (t.L1 - fc) + t.L2);
        }
        t.X0_cm = t.X0 / t.rho;
    }
    if (losses) {   // init.f:44-55 (trip_thru_target typeflag 4)
        TargetInitArgs a;
        a.using_Eloss = sw.using_Eloss;
        a.using_Coulomb = sw.using_Coulomb;
        a.the_cent = s.spec.e.theta;
        a.thp_cent = s.spec.p.theta;
        a.Pp_cent = s.spec.p.P;
        a.Mh2 = s.Mh2;
        a.Ebeam = s.Ebeam;
        a.Pe_cent = s.spec.e.P;
        losses->mostProbableLosses(a, t);
    }
    if (!sw.using_Eloss)
        for (auto& e : t.Eloss) e.ave = 0.0;   // init.f:56-60
    if (sw.using_Coulomb) {                    // init.f:62-82 (Aste et al. version)
        t.Coulomb.ave = 0.75 * 1.5 * (t.Z - 1.) * K::alpha * K::hbarc /
                        (1.1 * std::pow(t.A, 1. / 3.) + 0.86 * std::pow(t.A, -1. / 3.));
        t.Coulomb_constant = t.Coulomb.ave;
        t.Coulomb.min = t.Coulomb_constant;
        t.Coulomb.max = t.Coulomb_constant;
    } else {
        t.Coulomb.ave = 0.0;
        t.Coulomb_constant = 0.0;
        t.Coulomb.min = 0.0;
        t.Coulomb.max = 0.0;
    }

    // ---- dbase.f:555-558 theory_init
    if (f.doing_deuterium || (f.doing_heavy && !p.use_benhar_sf)) {
        std::string file;
        if (A == 2) file = "h2.theory";
        else if (A == 12) file = "c12.theory";
        else if (A == 56) file = "fe56.theory";
        else if (A == 197) file = "au197.theory";
        else {
            msg("No Theory File (spectral function) for A = " + std::to_string(t.A) + ", defaulting to c12.theory");
            file = "c12.theory";
        }
        const bool okT = ReadTheoryFile(joinPath(opt.data_dir, file), f.doing_heavy, s.theory, &s.messages);
        s.theory.file = file;
        if (!okT) {
            stop("THEORY_INIT failed!");
            return s;
        }
        // init.f:891 FORTRAN: theory_init overwrites doing_deuterium from the file contents
        s.flags.doing_deuterium = s.theory.nrhoPm == 1 && s.theory.E_Fermi < 1.0;
    }

    // ---- dbase.f:562-587 momentum distribution
    if (f.doing_deutpi || f.doing_hepi || f.doing_deutkaon || f.doing_hekaon || f.doing_deutsemi) {
        std::string file;
        if (f.doing_deutpi || f.doing_deutkaon || f.doing_deutsemi) file = "deut.dat";
        if (f.doing_hepi || f.doing_hekaon) {
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
            return s;
        }
        s.pdist.file = file;
    }

    // ---- dbase.f:594-621 spectral function (sf_lookup_init): not ported in P2.1
    if (f.doing_hepi || f.doing_hekaon || (f.doing_heavy && p.use_benhar_sf))
        s.not_ported.push_back("sf_lookup_init (dbase.f:594-621): Benhar spectral function; Emval(numEm) for the "
                               "VERTEXedge%Em%max of A>2 pion/kaon/delta/rho is left at its min_max_init value");
    // ---- dbase.f:624-681 kaon model tables: they belong to the kaon cross-section model
    if (f.doing_kaon) s.not_ported.push_back("kaon model tables weightc/weightd/saghai (dbase.f:624-681)");

    // =========================================================================================
    // init.f:limits_init (l.91-575), called at dbase.f:685
    // =========================================================================================
    Slops& slop = s.slop;
    struct SlopParam { double d, t, p; };
    // simulate.inc:17-36
    auto slopFor = [](int arm_code, SlopParam& sp) {
        switch (arm_code) {
            case 1: sp = {0.5, 0.005, 0.005}; return true;          // HMS
            case 2: sp = {1.0, 0.008, 0.008}; return true;          // SOS
            case 3: sp = {0.5, 0.005, 0.005}; return true;          // HRSR
            case 4: sp = {0.5, 0.005, 0.005}; return true;          // HRSL
            case 5: case 6: sp = {0.5, 0.005, 0.005}; return true;  // SHMS
            default: return false;                                  // 7/8 (calo): slop stays 0
        }
    };
    SlopParam sp;
    if (sw.using_E_arm_montecarlo && slopFor(p.electron_arm, sp)) {   // init.f:151-173
        slop.MC.e.delta.used = sp.d;
        slop.MC.e.yptar.used = sp.t;
        slop.MC.e.xptar.used = sp.p;
    }
    if (sw.using_P_arm_montecarlo && slopFor(p.hadron_arm, sp)) {     // init.f:174-196
        slop.MC.p.delta.used = sp.d;
        slop.MC.p.yptar.used = sp.t;
        slop.MC.p.xptar.used = sp.p;
    }

    // init.f:201-212 add slop to SPedges
    auto widen = [](ArmCuts& c, const Slops::Arm& a) {
        c.delta.min = c.delta.min - a.delta.used;
        c.delta.max = c.delta.max + a.delta.used;
        c.yptar.min = c.yptar.min - a.yptar.used;
        c.yptar.max = c.yptar.max + a.yptar.used;
        c.xptar.min = c.xptar.min - a.xptar.used;
        c.xptar.max = c.xptar.max + a.xptar.used;
    };
    widen(s.SPedge.e, slop.MC.e);
    widen(s.SPedge.p, slop.MC.p);

    // init.f:216-225 TRUE edges
    EdgeTrue& edge = s.edge;
    const double dE = p.dE_edge_test;
    const SpecArm& se = s.spec.e;
    const SpecArm& spp = s.spec.p;
    edge.e.E.min = (1. + s.SPedge.e.delta.min / 100.) * se.P + t.Coulomb.min - dE;
    edge.e.E.max = (1. + s.SPedge.e.delta.max / 100.) * se.P + t.Coulomb.max + dE;
    Cut pp;
    pp.min = (1. + s.SPedge.p.delta.min / 100.) * spp.P - dE;
    pp.max = (1. + s.SPedge.p.delta.max / 100.) * spp.P + dE;
    pp.min = std::max(0.001e0, pp.min);
    edge.p.E.min = std::sqrt(pp.min * pp.min + s.Mh2);
    edge.p.E.max = std::sqrt(pp.max * pp.max + s.Mh2);

    // init.f:228-246 extreme theta/z for extreme_trip_thru_target
    ExtremeTripArgs xa;
    xa.Ebeam = s.Ebeam;
    {
        const ArmCuts& a = s.SPedge.e;
        xa.the.max = std::acos((se.cos_th - se.sin_th * a.yptar.max) /
                               std::sqrt(1. + a.yptar.max * a.yptar.max + a.xptar.max * a.xptar.max));
        xa.the.min = std::acos((se.cos_th - se.sin_th * a.yptar.min) / std::sqrt(1. + a.yptar.min * a.yptar.min));
    }
    {
        const ArmCuts& a = s.SPedge.p;
        xa.thp.max = std::acos((spp.cos_th - spp.sin_th * a.yptar.max) /
                               std::sqrt(1. + a.yptar.max * a.yptar.max + a.xptar.max * a.xptar.max));
        xa.thp.min = std::acos((spp.cos_th - spp.sin_th * a.yptar.min) / std::sqrt(1. + a.yptar.min * a.yptar.min));
    }
    xa.z.min = -0.5 * t.length;
    xa.z.max = 0.5 * t.length;
    xa.pe = edge.e.E;
    xa.pp = pp;
    xa.m = s.Mh;
    if (losses) {
        losses->extremeLosses(xa, t);
        s.target_losses_from_provider = true;
    }
    if (!sw.using_Eloss) {   // init.f:247-252
        for (auto& e : t.Eloss) {
            e.min = 0.0;
            e.max = 0.0;
        }
    }
    if (!sw.mc_smear) t.musc_max = {0., 0., 0.};   // init.f:254-258

    // init.f:260-271
    edge.e.E.min = edge.e.E.min + t.Eloss[1].min;
    edge.e.E.max = edge.e.E.max + t.Eloss[1].max;
    edge.p.E.min = edge.p.E.min + t.Eloss[2].min;
    edge.p.E.max = edge.p.E.max + t.Eloss[2].max;
    edge.e.yptar.min = s.SPedge.e.yptar.min - t.musc_max[1];
    edge.e.yptar.max = s.SPedge.e.yptar.max + t.musc_max[1];
    edge.e.xptar.min = s.SPedge.e.xptar.min - t.musc_max[1];
    edge.e.xptar.max = s.SPedge.e.xptar.max + t.musc_max[1];
    edge.p.yptar.min = s.SPedge.p.yptar.min - t.musc_max[2];
    edge.p.yptar.max = s.SPedge.p.yptar.max + t.musc_max[2];
    edge.p.xptar.min = s.SPedge.p.xptar.min - t.musc_max[2];
    edge.p.xptar.max = s.SPedge.p.xptar.max + t.musc_max[2];

    // init.f:283 measured beam energy used in recon
    s.Ebeam_vertex_ave = s.Ebeam + t.Coulomb.ave - t.Eloss[0].ave;

    // init.f:287-298 slops in energies and momenta
    const double Ebeam_max = s.Ebeam + s.dEbeam / 2. - t.Eloss[0].min + t.Coulomb.max;
    const double Ebeam_min = s.Ebeam - s.dEbeam / 2. - t.Eloss[0].max + t.Coulomb.min;
    const double slop_Coulomb = t.Coulomb.max - t.Coulomb.ave;
    double slop_Ebeam = s.dEbeam / 2. + slop_Coulomb;
    double slop_Ee = slop.MC.e.delta.used / 100. * se.P + slop_Coulomb;
    const double r = std::sqrt(edge.p.E.max * edge.p.E.max - s.Mh2);
    double slop_Ep = std::sqrt((r + slop.MC.p.delta.used / 100. * spp.P) * (r + slop.MC.p.delta.used / 100. * spp.P) +
                               s.Mh2) -
                     edge.p.E.max;
    slop_Ebeam = slop_Ebeam + (t.Eloss[0].max - t.Eloss[0].min);
    slop_Ee = slop_Ee + (t.Eloss[1].max - t.Eloss[1].min);
    slop_Ep = slop_Ep + (t.Eloss[2].max - t.Eloss[2].min);

    if (f.doing_heavy) {   // init.f:300-305 'reconstructed' Em cuts
        slop.total.Em.used = slop_Ebeam + slop_Ee + slop_Ep + dE;
        edge.Em.min = s.cuts.Em.min - slop.total.Em.used;
        edge.Em.max = s.cuts.Em.max + slop.total.Em.used;
        edge.Em.min = std::max(0.e0, edge.Em.min);
    }

    // init.f:325-371 VERTEXedge Em/Pm
    EdgeTrue& V = s.VERTEXedge;
    const TheoryData& th = s.theory;
    auto pdistMax = [&]() { return s.pdist.pval.empty() ? 0.0 : s.pdist.pval.back(); };   // pval(nump)
    if (f.doing_hyd_elast) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    } else if (s.flags.doing_deuterium) {
        V.Em.min = K::Mp + K::Mn - t.M;
        V.Em.max = K::Mp + K::Mn - t.M;
        V.Pm.min = 0.0;
        V.Pm.max = th.Pm.empty() ? 0.0 : std::max(std::abs(th.Pm[0].min), std::abs(th.Pm[0].max));
    } else if (f.doing_heavy) {
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
        if (p.use_benhar_sf) {
            V.Pm.max = 790.0;
        } else {
            for (int i = 0; i < th.nrhoPm; ++i) {
                const double t1 = std::max(std::abs(th.Pm[i].min), std::abs(th.Pm[i].max));
                V.Pm.max = std::max(V.Pm.max, t1);
            }
        }
        V.Em.min = th.E_Fermi;
        V.Em.max = 1000.;   // Need Egamma_tot_max for good limit.
    } else if (f.doing_hydpi || f.doing_hydkaon || f.doing_hyddelta || f.doing_hydrho) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    } else if (f.doing_deutpi || f.doing_deutkaon || f.doing_deutdelta || f.doing_deutrho) {
        V.Em.min = K::Mp + K::Mn - t.M;
        V.Em.max = K::Mp + K::Mn - t.M;
        V.Pm.min = 0.0;
        V.Pm.max = pdistMax();
    } else if (f.doing_hepi || f.doing_hekaon || f.doing_hedelta || f.doing_herho) {
        V.Em.min = t.Mtar_struck + t.Mrec - t.M;
        // V.Em.max = Emval(numEm): sf_lookup_init not ported (see not_ported); keeps +1e10
        V.Pm.min = 0.0;
        V.Pm.max = pdistMax();
    } else if (f.doing_semi) {
        V.Em.min = 0.0;
        V.Em.max = 0.0;
        V.Pm.min = 0.0;
        V.Pm.max = 0.0;
    }
    msg("E_bind = " + std::to_string(V.Em.min) + " MeV in limits_init (QF only)");

    // init.f:379-392 recoiling (A-1) system
    if (f.doing_hyd_elast || f.doing_hydpi || f.doing_hydkaon || f.doing_hyddelta || f.doing_hydrho || f.doing_semi) {
        V.Mrec.min = 0.0;
        V.Mrec.max = 0.0;
        V.Trec.min = 0.0;
        V.Trec.max = 0.0;
    } else {
        V.Mrec.min = t.M - t.Mtar_struck + V.Em.min;
        V.Mrec.max = t.M - t.Mtar_struck + V.Em.max;
        V.Trec.min = std::sqrt(V.Mrec.max * V.Mrec.max + V.Pm.min * V.Pm.min) - V.Mrec.max;
        V.Trec.max = std::sqrt(V.Mrec.min * V.Mrec.min + V.Pm.max * V.Pm.max) - V.Mrec.min;
    }

    // init.f:399-406
    if (f.doing_eep || f.doing_semi) {
        V.Trec_struck.min = 0.;
        V.Trec_struck.max = 0.;
    } else {
        V.Trec_struck.min = 0.;
        V.Trec_struck.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.e.E.min - edge.p.E.min - V.Em.min -
                            V.Trec.min;
    }

    // init.f:415-432 radiation limits
    RadSetup& rad = s.rad;
    rad.Egamma_tot_max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.e.E.min - edge.p.E.min - V.Em.min -
                         V.Trec.min - V.Trec_struck.min;
    if (f.doing_heavy) {
        const double t2 = (edge.Em.max - V.Em.min) + (V.Trec.max - V.Trec.min);
        rad.Egamma_tot_max = std::min(rad.Egamma_tot_max, t2);
    }
    if (sw.hardwired_rad) rad.Egamma_tot_max = p.Egamma_gen_max;
    if (!sw.using_rad) rad.Egamma_tot_max = 0.0;
    if (sw.doing_tail[0]) rad.Egamma1_max = rad.Egamma_tot_max;
    if (sw.doing_tail[1]) rad.Egamma2_max = rad.Egamma_tot_max;
    if (sw.doing_tail[2]) rad.Egamma3_max = rad.Egamma_tot_max;
    if (f.doing_heavy) {
        V.Em.min = std::max(V.Em.min, edge.Em.min - rad.Egamma_tot_max);
        V.Em.max = std::min(V.Em.max, edge.Em.max);
    }

    // init.f:441-475 generated summed energies
    if (f.doing_hyd_elast) {
        gen.sumEgen.min = 0.0;
        gen.sumEgen.max = 0.0;
    } else if (f.doing_heavy) {
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - V.Trec.min - V.Em.min;
        gen.sumEgen.min = Ebeam_min + t.Mtar_struck - V.Trec.max - V.Em.max - rad.Egamma1_max;
        gen.sumEgen.max = std::min(gen.sumEgen.max, edge.e.E.max + edge.p.E.max + rad.Egamma_tot_max);
        gen.sumEgen.min = std::max(gen.sumEgen.min, edge.e.E.min + edge.p.E.min);
    } else if (f.doing_semi) {
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck;
        gen.sumEgen.min = edge.e.E.min + edge.p.E.min;
    } else {   // generated ELECTRON energy limits
        gen.sumEgen.max = Ebeam_max + t.Mtar_struck - t.Mrec_struck - edge.p.E.min - V.Em.min - V.Trec.min -
                          V.Trec_struck.min;
        gen.sumEgen.min = Ebeam_min + t.Mtar_struck - t.Mrec_struck - edge.p.E.max - V.Em.max - V.Trec.max -
                          V.Trec_struck.max - rad.Egamma_tot_max;
        gen.sumEgen.max = std::min(gen.sumEgen.max, edge.e.E.max + rad.Egamma2_max);
        gen.sumEgen.min = std::max(gen.sumEgen.min, edge.e.E.min);
    }
    gen.sumEgen.min = gen.sumEgen.min - dE;
    gen.sumEgen.max = gen.sumEgen.max + dE;
    gen.sumEgen.min = std::max(0.e0, gen.sumEgen.min);

    // init.f:480-498 E arm generation limits
    if (f.doing_hyd_elast) {
        gen.e.E.min = edge.e.E.min;
        gen.e.E.max = edge.e.E.max + rad.Egamma2_max;
    } else if (s.flags.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_rho || f.doing_delta) {
        gen.e.E.min = gen.sumEgen.min;
        gen.e.E.max = gen.sumEgen.max;
    } else if (f.doing_heavy || f.doing_semi) {
        gen.e.E.min = gen.sumEgen.min - edge.p.E.max - rad.Egamma3_max;
        gen.e.E.max = gen.sumEgen.max - edge.p.E.min;
    }
    gen.e.E.min = std::max(gen.e.E.min, edge.e.E.min);
    gen.e.E.max = std::min(gen.e.E.max, edge.e.E.max + rad.Egamma2_max);
    gen.e.delta.min = (gen.e.E.min / se.P - 1.) * 100.;
    gen.e.delta.max = (gen.e.E.max / se.P - 1.) * 100.;
    gen.e.yptar.min = edge.e.yptar.min;
    gen.e.yptar.max = edge.e.yptar.max;
    gen.e.xptar.min = edge.e.xptar.min;
    gen.e.xptar.max = edge.e.xptar.max;

    // init.f:503-517 P arm generation limits
    if (f.doing_hyd_elast || s.flags.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_rho || f.doing_delta) {
        gen.p.E.min = edge.p.E.min;
        gen.p.E.max = edge.p.E.max + rad.Egamma3_max;
    } else if (f.doing_heavy || f.doing_semi) {
        gen.p.E.min = gen.sumEgen.min - edge.e.E.max - rad.Egamma2_max;
        gen.p.E.max = gen.sumEgen.max - edge.e.E.min;
    }
    gen.p.E.min = std::max(gen.p.E.min, edge.p.E.min);
    gen.p.E.max = std::min(gen.p.E.max, edge.p.E.max + rad.Egamma3_max);
    gen.p.delta.min = (std::sqrt(gen.p.E.min * gen.p.E.min - s.Mh2) / spp.P - 1.) * 100.;
    gen.p.delta.max = (std::sqrt(gen.p.E.max * gen.p.E.max - s.Mh2) / spp.P - 1.) * 100.;
    gen.p.yptar.min = edge.p.yptar.min;
    gen.p.yptar.max = edge.p.yptar.max;
    gen.p.xptar.min = edge.p.xptar.min;
    gen.p.xptar.max = edge.p.xptar.max;

    // init.f:520-535 axis specs for the diagnostic histograms (H%gen = H%geni = H%RECON)
    constexpr double nHbins = 50;   // histograms.inc
    DiagHistAxes& H = s.hist;
    H.e.delta.min = gen.e.delta.min;
    H.e.yptar.min = gen.e.yptar.min;
    H.e.xptar.min = -gen.e.xptar.max;
    H.e.delta.bin = (gen.e.delta.max - gen.e.delta.min) / nHbins;
    H.e.yptar.bin = (gen.e.yptar.max - gen.e.yptar.min) / nHbins;
    H.e.xptar.bin = (gen.e.xptar.max - gen.e.xptar.min) / nHbins;
    H.p.delta.min = gen.p.delta.min;
    H.p.yptar.min = gen.p.yptar.min;
    H.p.xptar.min = -gen.p.xptar.max;
    H.p.delta.bin = (gen.p.delta.max - gen.p.delta.min) / nHbins;
    H.p.yptar.bin = (gen.p.yptar.max - gen.p.yptar.min) / nHbins;
    H.p.xptar.bin = (gen.p.xptar.max - gen.p.xptar.min) / nHbins;
    H.Em.min = V.Em.min;
    H.Em.bin = (std::max(100.e0, V.Em.max) - V.Em.min) / nHbins;
    H.Pm.min = V.Pm.min;
    H.Pm.bin = (std::max(100.e0, V.Pm.max) - V.Pm.min) / nHbins;

    // =========================================================================================
    // dbase.f:689-958 announcements: only the STOP branches matter for the setup
    // =========================================================================================
    const ReactionFlags& ff = s.flags;
    auto ngammaOk = [&]() { return p.pizero_ngamma == 1 || p.pizero_ngamma == 2; };
    if (ff.doing_eep) {
        if (!(ff.doing_hyd_elast || ff.doing_deuterium || ff.doing_heavy)) {
            stop("I don't have ANY idea what (e,e'p) we're doing!!!");
            return s;
        }
    } else if (ff.doing_semi) {
        if (ff.doing_semipi) {
            if ((ff.doing_hydsemi || ff.doing_deutsemi || ff.doing_hesemi) && ff.doing_pizero && !ngammaOk()) {
                stop("pizero_ngamma not set correctly(must be 1 or 2), stopping");
                return s;
            }
        } else if (!ff.doing_semika) {
            stop("I don't have ANY idea what (e,e'm)X we're doing!!!");
            return s;
        }
    } else if (ff.doing_delta) {
        // FORTRAN-QUIRK: dbase.f:780-788 tests doing_hyddelta, doing_deutpi, doing_hepi (not doing_deutdelta /
        // doing_hedelta), so D(e,e'p)pi and A(e,e'p)pi runs STOP here.
        if (!(ff.doing_hyddelta || ff.doing_deutpi || ff.doing_hepi)) {
            stop("I don't have ANY idea what (e,e'p)pi we're doing!!!");
            return s;
        }
    } else if (ff.doing_pion) {
        if (!(ff.doing_hydpi || ff.doing_deutpi || ff.doing_hepi)) {
            stop("I don't have ANY idea what (e,e'pi) we're doing!!!");
            return s;
        }
        if (ff.doing_pizero && !ngammaOk()) {
            stop("pizero_ngamma not set correctly(must be 1 or 2), stopping");
            return s;
        }
        // dbase.f:823-826. FORTRAN-QUIRK: this also stops charged H(e,e'pi-)Delta++ (which_pion=3 without
        // doing_pizero), although the message only talks about pi0 from a neutron.
        if (p.which_pion == 3 && ff.doing_hydpi) {
            stop("Can't do Delta+ production from neutron for pi0");
            return s;
        }
    } else if (ff.doing_kaon) {
        if (!(ff.doing_hydkaon || ff.doing_deutkaon || ff.doing_hekaon)) {
            stop("I don't have ANY idea what (e,e'K) we're doing!!!");
            return s;
        }
        const int w = p.which_kaon;
        if (!(w == 0 || w == 1 || w == 2 || w == 10 || w == 11 || w == 12)) {
            stop("I don't have ANY idea what (e,e'K) we're doing!!!");
            return s;
        }
    } else if (ff.doing_rho) {
        if (ff.doing_hydrho) {
            if (p.targ.A >= 3) {
                stop("A(e,e'rho) Not yet implemented");
                return s;
            }
        } else if (!(ff.doing_deutrho || ff.doing_herho)) {
            stop("I don't have ANY idea what (e,e'rho) we're doing!!!");
            return s;
        }
    } else if (!ff.doing_phsp) {
        stop("I don't have ANY idea what we're doing!!!");
        return s;
    }
    if (p.electron_arm == 6 || p.hadron_arm == 6) {   // dbase.f:895-897, 915-917
        stop("This spectrometer option no longer available");
        return s;
    }
    if (ff.doing_phsp) msg("**WARNING: phsp option is not really implemented - buyer beware!!");

    // =========================================================================================
    // init.f:radc_init (l.576-654), called at simc.f:84
    // =========================================================================================
    if (sw.extrad_flag == 0) {
        if (sw.rad_flag == 0) {
            sw.extrad_flag = 3;
        } else if (sw.rad_flag == 1 || sw.rad_flag == 2 || sw.rad_flag == 3) {
            sw.extrad_flag = 1;
        }
    } else if (sw.extrad_flag < 0) {
        stop("Imbecile! check your stupid setting of EXTRAD_FLAG");
        return s;
    }
    rad.etatzai = (12.0 + (t.Z + 1.) / (t.Z * t.L1 + t.L2)) / 9.0;
    sw.produce_output = p.debug[0];
    sw.exponentiate = p.use_expon == 1;
    sw.include_hard = true;
    sw.calculate_spence = true;

    // =========================================================================================
    // simc.f:94-101 target factor and luminosity
    // =========================================================================================
    s.targetfac = t.mass_amu / 3.75914e+6 / (t.abundancy / 100.) * std::abs(std::cos(t.angle)) / (t.thick * 1000.);
    s.luminosity = p.EXPER.charge / s.targetfac;
    return s;
}

}  // namespace simc
