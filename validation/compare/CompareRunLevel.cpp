// validation/compare/CompareRunLevel.cpp — the run-level tests of simc_compare: pulls of the .hist/summary numbers
// (counters, wtcontr, averages, normalisation) and of the .geni stop counters, with their statistical errors.
#include "CompareCommon.h"

namespace simccmp {

void RunLevelTests(Session& st) {
    Config& cfg = st.cfg;
    std::vector<std::string>& warnings = st.warnings;
    RunSummary &RA = st.RA, &RB = st.RB;
    TreeData &TA = st.TA, &TB = st.TB;
    Norm &NA = st.NA, &NB = st.NB;
    const std::vector<double>* WA = st.WA;
    const std::vector<double>* WB = st.WB;
    const bool haveCounters = st.haveCounters;
    auto add = [&](Test t) { return st.add(std::move(t)); };

    if (haveCounters) {
        if (RA.ncontribute != int64_t(TA.n)) warnings.push_back("A: tree entries " + std::to_string(TA.n) + " != Ncontribute " + std::to_string(RA.ncontribute));
        if (RB.ncontribute != int64_t(TB.n)) warnings.push_back("B: tree entries " + std::to_string(TB.n) + " != Ncontribute " + std::to_string(RB.ncontribute));
        Test t;
        t.category = "norm";
        t.name = "acceptance Ncontribute/Ntried";
        t.method = "binomial";
        binomialPull(RA.ncontribute, RA.ntried, RB.ncontribute, RB.ntried, t.a, t.b, t.ea, t.eb, t.pull);
        t.p = pullPValue(t.pull);
        add(t);
        Test u;
        u.category = "norm";
        u.name = "Nco_no_rad_prot/Ncontribute";
        u.method = "binomial";
        if (binomialPull(RA.nco_no_rad_prot, RA.ncontribute, RB.nco_no_rad_prot, RB.ncontribute, u.a, u.b, u.ea, u.eb, u.pull)) {
            u.p = pullPValue(u.pull);
            add(u);
        }
    } else {
        warnings.push_back("run-level counters missing in A or B: normalisation pulls skipped");
    }

    // Tree yield per tried event: Y = scale*sum W ; see README
    auto treeYield = [](const std::vector<double>& W, const Norm& n, double& Y, double& sY) {
        double S1 = 0, S2 = 0;
        for (double w : W) if (std::isfinite(w)) { S1 += w; S2 += w * w; }
        if (n.perTried) {
            perTriedMean(S1, S2, n.ntried, Y, sY);
            double f = n.scale * n.ntried;  // 1 for per_tried, lumi*genvol for normfac
            Y *= f;
            sY *= f;
        } else {
            Y = n.scale * S1;
            sY = n.scale * std::sqrt(S2);
        }
    };
    double YA, sYA, YB, sYB;
    treeYield(*WA, NA, YA, sYA);
    treeYield(*WB, NB, YB, sYB);
    {
        Test t;
        t.category = "norm";
        t.name = "yield sum(" + RA.weightBranch + ") [" + NA.desc + "]";
        t.method = "pull";
        t.a = YA; t.ea = sYA; t.b = YB; t.eb = sYB;
        t.pull = pull(YA, sYA, YB, sYB);
        t.p = pullPValue(t.pull);
        add(t);
    }
    // wtcontr (.hist): sum of main%weight over the pass_cuts events only (simc.f:229-241, 311-313), times normfac/nevent.
    // D-053: its relative error is that of the pass_cuts subset of the tree, sqrt(sum w^2)/sum w of the events in
    // the run's own pass_cuts window (perTriedMean over Ntried as the tree yield); without the window data (SPedge, slop,
    // electron_arm, the mode) the relative error of the whole tree yield is used, with a note.
    if (const Scalar *wa = RA.scalar("wtcontr"), *wb = RB.scalar("wtcontr"); wa && wb) {
        Test t;
        t.category = "norm";
        t.name = "wtcontr";
        t.method = "pull";
        t.a = wa->value; t.b = wb->value;
        std::string why[2];
        auto relerr = [&](int k) -> double {
            const RunSummary& R = k ? RB : RA;
            const TreeData& T = k ? TB : TA;
            const Norm& N = k ? NB : NA;
            const std::vector<double>& W = k ? *WB : *WA;
            const double Yall = k ? YB : YA, sYall = k ? sYB : sYA;
            auto fallback = [&](const std::string& reason) {
                why[k] = reason;
                return Yall > 0 ? sYall / Yall : 0.0;
            };
            if (!R.has_spedge || !R.has_slop) return fallback("no SPedge/slop");
            if (R.fortran_compat < 0) return fallback("pass_cuts mode unknown (--mode)");
            const bool right = R.electron_arm == 1 || R.electron_arm == 3 || R.electron_arm == 7;
            const bool left = R.electron_arm == 2 || R.electron_arm == 4 || R.electron_arm == 5 || R.electron_arm == 8;
            if (!right && !left) return fallback("electron_arm unknown");
            const char* ep[2] = {right ? "hs" : "ss", right ? "ss" : "hs"};   // e, p prefixes
            const char* var[3] = {"delta", "yptar", "xptar"};
            const std::vector<double>* c[6];
            for (int a = 0; a < 2; ++a)
                for (int v = 0; v < 3; ++v) {
                    c[3 * a + v] = T.col(std::string(ep[a]) + var[v]);
                    if (!c[3 * a + v]) return fallback(std::string("no branch ") + ep[a] + var[v]);
                }
            // the window: printed SPedge (already + slop) -/+ slop used; angles: slop in mr
            double lo[6], hi[6];
            for (int i = 0; i < 6; ++i) {
                const double sl = (i % 3 == 0) ? R.slop_used[i] : R.slop_used[i] / 1000.;
                lo[i] = R.spedge_min[i] + sl;
                hi[i] = R.spedge_max[i] - sl;
            }
            if (R.fortran_compat == 1) hi[3] = R.spedge_max[0] - R.slop_used[3];   // simc.f:237 (D-006)
            double S1 = 0, S2 = 0;
            for (size_t i = 0; i < W.size(); ++i) {
                bool pass = std::isfinite(W[i]);
                for (int k6 = 0; k6 < 6 && pass; ++k6) {
                    const double x = (*c[k6])[i];
                    if (x <= lo[k6] || x >= hi[k6]) pass = false;
                }
                if (pass) { S1 += W[i]; S2 += W[i] * W[i]; }
            }
            if (!(S1 > 0)) return fallback("empty pass_cuts subset");
            double Y = 0, sY = 0;
            if (N.perTried) perTriedMean(S1, S2, N.ntried, Y, sY);
            else { Y = S1; sY = std::sqrt(S2); }
            // consistency: S1*normfac/nevent should be wtcontr
            if (const Scalar* nf = R.scalar("normfac"); nf && R.has_counters && R.nevent() > 0) {
                const double pred = S1 * nf->value / double(R.nevent());
                const double wt = (k ? wb : wa)->value;
                if (std::fabs(pred - wt) > 2e-3 * std::fabs(wt))
                    warnings.push_back(std::string(k ? "B" : "A") + ": pass_cuts subset of the tree gives wtcontr " +
                                       fmt(pred) + " != " + fmt(wt) + " (window mode " +
                                       (R.fortran_compat ? "compat" : "corrected") + ")");
            }
            why[k] = std::string("pass_cuts subset, ") + (R.fortran_compat ? "compat" : "corrected") + " window";
            return Y > 0 ? sY / Y : 0.0;
        };
        const double ra = relerr(0), rb = relerr(1);
        t.ea = std::fabs(wa->value) * ra;
        t.eb = std::fabs(wb->value) * rb;
        t.pull = pull(t.a, t.ea, t.b, t.eb);
        t.p = pullPValue(t.pull);
        t.note = "sigma = wtcontr * relerr: A " + why[0] + "; B " + why[1];
        for (int k = 0; k < 2; ++k)
            if (why[k].rfind("pass_cuts subset", 0) != 0)
                warnings.push_back(std::string(k ? "B" : "A") + ": wtcontr sigma from the whole tree yield (" + why[k] + ")");
        add(t);
        // info: consistency of wtcontr with the tree (pass_cuts subset, simc.f:311-313)
        for (int k = 0; k < 2; ++k) {
            const RunSummary& R = k ? RB : RA;
            const std::vector<double>& W = k ? *WB : *WA;
            const Scalar* nf = R.scalar("normfac");
            if (!nf || !R.has_counters || R.nevent() == 0) continue;
            double S = 0;
            for (double w : W) if (std::isfinite(w)) S += w;
            double pred = S * nf->value / double(R.nevent());
            double wt = (k ? wb : wa)->value;
            if (std::fabs(pred - wt) > 1e-3 * std::fabs(wt))
                warnings.push_back(std::string(k ? "B" : "A") + ": sum(Weight)*normfac/nevent = " + fmt(pred) +
                                   " differs from wtcontr = " + fmt(wt) + " (tree holds all successes, wtcontr only pass_cuts)");
        }
    }
    // AVERAGE.sigcc = sum(sigcc over successes)/nevent (simc.f:248,959). The pion/kaon/delta ntuple calls the same
    // quantity "siglab" (results_write.f:146, NtupleInit.f:129): used when the configured branch is absent.
    std::string sigcc_b = cfg.sigccBranch;
    if ((!TA.col(sigcc_b) || !TB.col(sigcc_b)) && TA.col("siglab") && TB.col("siglab")) sigcc_b = "siglab";
    if (const Scalar *sa = RA.scalar("sigcc_average"), *sb = RB.scalar("sigcc_average");
        sa && sb && haveCounters && TA.col(sigcc_b) && TB.col(sigcc_b)) {
        auto err = [&](const RunSummary& R, const std::vector<double>& x) {
            double S1 = 0, S2 = 0, n = 0;
            for (double v : x) if (std::isfinite(v)) { S1 += v; S2 += v * v; n += 1; }
            if (R.ngen < 0) {  // nevent = Ntried
                double Y, sY;
                perTriedMean(S1, S2, double(R.ntried), Y, sY);
                return sY;
            }
            if (n < 2) return 0.0;
            double m = S1 / n;
            return std::sqrt(std::max(0.0, S2 / n - m * m) / n);
        };
        Test t;
        t.category = "norm";
        t.name = "AVERAGE.sigcc";
        t.method = "pull";
        t.a = sa->value; t.b = sb->value;
        t.ea = std::hypot(err(RA, *TA.col(sigcc_b)), sa->halfulp);
        t.eb = std::hypot(err(RB, *TB.col(sigcc_b)), sb->halfulp);
        t.pull = pull(t.a, t.ea, t.b, t.eb);
        t.p = pullPValue(t.pull);
        add(t);
    }
    // normfac: random through Ntried when ngen>0 (info only: same content as acceptance)
    if (const Scalar *na = RA.scalar("normfac"), *nb = RB.scalar("normfac"); na && nb && haveCounters) {
        Test t;
        t.category = (RA.ngen > 0 && RB.ngen > 0) ? "info" : "det";
        t.name = "normfac";
        t.a = na->value; t.b = nb->value;
        if (t.category == "info") {
            t.method = "pull";
            auto rel = [](const RunSummary& R) {
                double a = double(R.ncontribute) / double(R.ntried);
                return std::sqrt(std::max(0.0, 1 - a) / double(R.ncontribute));
            };
            t.ea = t.a * rel(RA);
            t.eb = t.b * rel(RB);
            t.pull = pull(t.a, t.ea, t.b, t.eb);
            t.p = pullPValue(t.pull);
            t.counted = false;
            t.note = "Var(Ntried)/Ntried^2=(1-acc)/Ncontribute; not counted (duplicates acceptance)";
            add(t);
        }
    }
    // deterministic scalars
    {
        std::vector<std::string> det = {"sigcc_central", "charge_mC", "targetfac", "luminosity_ub",
                                        "luminosity_GeV2", "genvol"};
        if (haveCounters && !(RA.ngen > 0 && RB.ngen > 0)) det.push_back("normfac");
        for (const auto& k : det) {
            const Scalar *a = RA.scalar(k), *b = RB.scalar(k);
            if (!a || !b) continue;
            Test t;
            t.category = "det";
            t.name = k;
            t.method = "det";
            t.a = a->value; t.b = b->value;
            double tol = cfg.detRelTol * std::max(std::fabs(t.a), std::fabs(t.b)) + a->halfulp + b->halfulp;
            t.stat = std::fabs(t.a - t.b);
            t.ea = tol;
            t.pass = t.stat <= tol;
            t.p = t.pass ? 1.0 : 0.0;
            t.note = "tol " + fmt(tol, 3);
            add(t);
        }
        if (haveCounters) {
            Test t;
            t.category = "det";
            t.name = "Ngen (request)";
            t.method = "det";
            t.a = double(RA.ngen); t.b = double(RB.ngen);
            t.stat = std::fabs(t.a - t.b);
            t.pass = RA.ngen == RB.ngen;
            t.p = t.pass ? 1.0 : 0.0;
            add(t);
        }
    }
    // stop counters
    if (haveCounters) {
        static const std::set<std::string> nonCategory = {"trials", "hut", "q2_mid_dup"};
        for (const auto& [spec, mA] : RA.stops) {
            auto itB = RB.stops.find(spec);
            if (itB == RB.stops.end()) {
                warnings.push_back("stop counters for " + spec + " only in A");
                continue;
            }
            const auto& mB = itB->second;
            // Only the counters reported by BOTH runs are compared: the Fortran .geni does not print every
            // counter (SHMS q3_mid, D1_men, targ_*, cer1/2, BP_in/out; OUTPUT_FORMAT.md), while the SIMC++ summary
            // carries all of them; a counter missing on one side is not a measured 0. Such keys are listed as warnings
            // when they are non-zero; their counts end up in the 'unaccounted' cell of the partition test.
            std::vector<std::string> keys;
            const auto& oB = RB.stopOrder.at(spec);
            for (const auto& k : RA.stopOrder.at(spec))
                if (std::find(oB.begin(), oB.end(), k) != oB.end()) keys.push_back(k);
            for (int side = 0; side < 2; ++side) {
                const auto& order = side ? oB : RA.stopOrder.at(spec);
                const auto& mm = side ? itB->second : mA;
                for (const auto& k : order) {
                    if (std::find(keys.begin(), keys.end(), k) != keys.end()) continue;
                    auto it = mm.find(k);
                    if (it != mm.end() && it->second != 0)
                        warnings.push_back(spec + "." + k + " = " + std::to_string(it->second) + " only in " +
                                           (side ? "B" : "A") + " (not reported by the other run; not compared)");
                }
            }
            auto get = [](const std::map<std::string, int64_t>& m, const std::string& k) {
                auto it = m.find(k);
                return it == m.end() ? int64_t(0) : it->second;
            };
            if (cfg.stopPerCounter) {
                for (const auto& k : keys) {
                    Test t;
                    t.category = "stop";
                    t.name = spec + "." + k + "/Ntried";
                    t.method = "binomial";
                    if (!binomialPull(get(mA, k), RA.ntried, get(mB, k), RB.ntried, t.a, t.b, t.ea, t.eb, t.pull)) continue;
                    t.p = pullPValue(t.pull);
                    add(t);
                }
            }
            // homogeneity of the partition of the spectrometer trials
            std::vector<int64_t> ca, cb;
            int64_t sa = 0, sb = 0;
            for (const auto& k : keys) {
                if (nonCategory.count(k)) continue;
                ca.push_back(get(mA, k));
                cb.push_back(get(mB, k));
                sa += ca.back();
                sb += cb.back();
            }
            int64_t ua = get(mA, "trials") - sa, ub = get(mB, "trials") - sb;
            if (ua < 0 || ub < 0) {
                warnings.push_back(spec + ": stop counters exceed trials; partition test skipped");
                continue;
            }
            ca.push_back(ua);  // 'unaccounted' (e.g. SHMS q3_mid, never printed: simc.f:520)
            cb.push_back(ub);
            Test t;
            t.category = "stop";
            t.name = spec + " partition of trials";
            t.method = "contingency";
            t.a = double(get(mA, "trials"));
            t.b = double(get(mB, "trials"));
            if (contingencyChi2(ca, cb, t.stat, t.ndf, t.p)) add(t);
        }
        for (const auto& [spec, m] : RB.stops)
            if (!RA.stops.count(spec)) warnings.push_back("stop counters for " + spec + " only in B");
    }

}

}  // namespace simccmp
