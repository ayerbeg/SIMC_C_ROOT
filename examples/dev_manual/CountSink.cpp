// examples/dev_manual/CountSink.cpp — the IOutputSink + IParallelSink example of the developer manual (task P15.D1,
// docs/manual/developer/04_extending.md). "example_count": counts the tries and the successes and sums the weights of
// the successes; writes "<output_dir>/<run_name>.count.txt" at the end of the run. It works with --threads N: every
// worker fills its own CountSink, the run's sink adds them up in the worker order.
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/RunSummary.h"
#include "simc/stages/IOutputSink.h"
#include "simc/stages/IParallelSink.h"

namespace mylab {

class CountSink : public simc::IOutputSink, public simc::IParallelSink {
public:
    std::string name() const override { return "example_count"; }
    void configure(const simc::Json& config, const simc::SimcParameters&) override {
        if (config.is_object()) file_ = config.value("file", file_);
    }
    void init(const simc::StageContext&) override {}
    void beginRun(const simc::RunSetup& setup) override {
        if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".count.txt";
    }
    void fill(const simc::MainEvent& main, const simc::Event&, const simc::Event&, const simc::Event&,
              bool success) override {
        ++tries_;
        if (success) {
            ++successes_;
            sumWeight_ += main.weight;
        }
    }
    void endRun(const simc::RunSummary& summary) override {
        std::ofstream f(file_);
        if (!f) throw std::runtime_error("example_count: cannot write " + file_);
        f << "tries " << tries_ << "\nsuccesses " << successes_ << "\nsum_weight " << sumWeight_ << "\nntried "
          << summary.ntried << "\n";
    }

    // --threads N: a fresh sink per worker (fill only), merged in the worker order 0..N-1 before endRun().
    std::unique_ptr<simc::IOutputSink> makeWorkerSink(int, int) const override {
        auto w = std::make_unique<CountSink>();
        w->file_ = file_;   // configured like this one; a worker sink never writes (no endRun)
        return w;
    }
    void mergeWorkerSink(simc::IOutputSink& sink) override {
        const auto& w = dynamic_cast<const CountSink&>(sink);
        tries_ += w.tries_;
        successes_ += w.successes_;
        sumWeight_ += w.sumWeight_;
    }

    long tries() const { return tries_; }
    long successes() const { return successes_; }
    double sumWeight() const { return sumWeight_; }

private:
    std::string file_;
    long tries_{0}, successes_{0};
    double sumWeight_{0};
};

}  // namespace mylab

SIMC_REGISTER(simc::IOutputSink, "example_count", mylab::CountSink)
