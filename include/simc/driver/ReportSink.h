/**
 * @file
 * @ingroup driver
 * @brief simc::driver::IReportSink and simc::driver::FortranHistSink: sinks of the full end-of-run report (.hist, .geni, .gen).
 *
 * Optional extension of IOutputSink for sinks that need the full end-of-run report of the driver (RunReport: the
 * central event, contrib/slop ranges, resolutions ...), not only the RunSummary. SimcRun calls writeReport() after
 * endRun() on every sink that also implements IReportSink (task P7.2).
 *
 * FortranHistSink ("fortran_hist"): the three text files of simc.f in the exact gfortran layout:
 *   <run>.hist  simc.f:report (l.644-1142)
 *   <run>.geni  simc.f l.462-532 (stop counters)
 *   <run>.gen   simc.f l.534-612 (the H%RECON, H%gen, H%geni histograms, 50 bins, histograms.inc), filled in fill()
 *               exactly like simc.f:250-287 with inc() (simc.f:624-640). The .gen quirks are kept: the "Original E arm"
 *               block prints H%gen%e%delta but H%geni%e%yptar/xptar, all on the H%gen bin centres; H%gen%Pm is never
 *               filled; H%RECON xptar is not negated.
 * JSON config (optional): {"hist": "<path>", "geni": "<path>", "gen": "<path>"}.
 */
#ifndef SIMC_DRIVER_REPORTSINK_H
#define SIMC_DRIVER_REPORTSINK_H

#include <array>
#include <memory>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/driver/RunReport.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IParallelSink.h"

namespace simc {
namespace driver {

/// @addtogroup driver
/// @{

/// @brief Optional extension of IOutputSink: a sink that also receives the full driver report (RunReport).
///
/// SimcRun calls writeReport() after endRun() on every sink that implements this interface.
class IReportSink {
public:
    /// Virtual destructor (interface).
    virtual ~IReportSink() = default;
    /// @brief Write the end-of-run report.
    /// @param setup the run set-up.
    /// @param report the driver report (central event, totals, normalisation, stop counters).
    virtual void writeReport(const RunSetup& setup, const RunReport& report) = 0;
};

/// @brief Sink "fortran_hist": the .hist, .geni and .gen text files in the exact gfortran layout.
///
/// Not thread-safe; in a threaded run every worker has its own copy (makeWorkerSink) and the copies are merged.
class FortranHistSink : public IOutputSink, public IReportSink, public IParallelSink {
public:
    static constexpr int kNHbins = 50;   ///< histograms.inc nHbins

    /// Fortran: histograms.inc hist_entry (min, bin, buf(nHbins)); buf[1..50] used.
    struct Hist {
        HistAxis ax;  ///< histogram axis (min, bin width).
        std::array<double, kNHbins + 1> buf{};  ///< bin contents, 1-based (buf[0] unused).
    };
    /// Fortran: type histograms (H%RECON, H%gen, H%geni): e/p delta, yptar, xptar, Em, Pm.
    struct HistSet {
        /// Histograms of one arm.
        struct Arm { Hist delta, yptar, xptar; /**< delta, yptar, xptar histograms */ } e, p;  ///< histograms of one arm; e = electron arm, p = hadron arm.
        Hist Em, Pm;  ///< missing energy and missing momentum histograms.
    };

    std::string name() const override { return "fortran_hist"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext&) override {}
    void beginRun(const RunSetup& setup) override;
    void fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon, bool success) override;
    void endRun(const RunSummary&) override {}
    void writeReport(const RunSetup& setup, const RunReport& report) override;
    /// P15.4 (--threads N): worker histograms in memory; merge = the 50-bin buffers summed (worker order).
    /// @param worker worker index (0-based).
    /// @param nworkers number of workers.
    /// @return a new in-memory sink for that worker.
    std::unique_ptr<IOutputSink> makeWorkerSink(int worker, int nworkers) const override;
    void mergeWorkerSink(IOutputSink& worker) override;

    /// Fortran: simc.f:inc(hist, val, weight) (l.624-640).
    /// @param h the histogram.
    /// @param val value to bin.
    /// @param weight weight added to the bin.
    static void Inc(Hist& h, double val, double weight);
    /// The .gen text (simc.f l.534-612).
    /// @return the .gen file text.
    std::string formatGen() const;

    /// @return the H%RECON histograms.
    const HistSet& recon() const { return RECON_; }
    /// @return the H%gen histograms.
    const HistSet& gen() const { return gen_; }
    /// @return the H%geni histograms.
    const HistSet& geni() const { return geni_; }
    /// @return the path of the .hist file.
    const std::string& histFile() const { return hist_; }
    /// @return the path of the .geni file.
    const std::string& geniFile() const { return geniFile_; }
    /// @return the path of the .gen file.
    const std::string& genFile() const { return genFile_; }

private:
    bool rhoGenFromOrig_{false};   ///< D-057 F18
    std::string hist_, geniFile_, genFile_;
    HistSet RECON_, gen_, geni_;
};

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_REPORTSINK_H
