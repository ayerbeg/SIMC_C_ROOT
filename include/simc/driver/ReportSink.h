// include/simc/driver/ReportSink.h
//
// Optional extension of IOutputSink for sinks that need the full end-of-run report of the driver (RunReport: the
// central event, contrib/slop ranges, resolutions ...), not only the RunSummary. SimcRun calls writeReport() after
// endRun() on every sink that also implements IReportSink (task P7.2).
//
// FortranHistSink ("fortran_hist"): the three text files of simc.f in the exact gfortran layout:
//   <run>.hist  simc.f:report (l.644-1142)
//   <run>.geni  simc.f l.462-532 (stop counters)
//   <run>.gen   simc.f l.534-612 (the H%RECON, H%gen, H%geni histograms, 50 bins, histograms.inc), filled in fill()
//               exactly like simc.f:250-287 with inc() (simc.f:624-640). The .gen quirks are kept: the "Original E arm"
//               block prints H%gen%e%delta but H%geni%e%yptar/xptar, all on the H%gen bin centres; H%gen%Pm is never
//               filled; H%RECON xptar is not negated.
// JSON config (optional): {"hist": "<path>", "geni": "<path>", "gen": "<path>"}.
#ifndef SIMC_DRIVER_REPORTSINK_H
#define SIMC_DRIVER_REPORTSINK_H

#include <array>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/driver/RunReport.h"
#include "simc/stages/IOutputSink.h"

namespace simc {
namespace driver {

class IReportSink {
public:
    virtual ~IReportSink() = default;
    virtual void writeReport(const RunSetup& setup, const RunReport& report) = 0;
};

class FortranHistSink : public IOutputSink, public IReportSink {
public:
    static constexpr int kNHbins = 50;   ///< histograms.inc nHbins

    /// Fortran: histograms.inc hist_entry (min, bin, buf(nHbins)); buf[1..50] used.
    struct Hist {
        HistAxis ax;
        std::array<double, kNHbins + 1> buf{};
    };
    /// Fortran: type histograms (H%RECON, H%gen, H%geni): e/p delta, yptar, xptar, Em, Pm.
    struct HistSet {
        struct Arm { Hist delta, yptar, xptar; } e, p;
        Hist Em, Pm;
    };

    std::string name() const override { return "fortran_hist"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(RunSetup&) override {}
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon, bool success) override;
    void endRun(const RunSummary&) override {}
    void writeReport(const RunSetup& setup, const RunReport& report) override;

    /// Fortran: simc.f:inc(hist, val, weight) (l.624-640).
    static void Inc(Hist& h, double val, double weight);
    /// The .gen text (simc.f l.534-612).
    std::string formatGen() const;

    const HistSet& recon() const { return RECON_; }
    const HistSet& gen() const { return gen_; }
    const HistSet& geni() const { return geni_; }
    const std::string& histFile() const { return hist_; }
    const std::string& geniFile() const { return geniFile_; }
    const std::string& genFile() const { return genFile_; }

private:
    std::string hist_, geniFile_, genFile_;
    HistSet RECON_, gen_, geni_;
};

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_REPORTSINK_H
