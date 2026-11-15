// src/output/RootTreeSink.cpp
// Fortran: results_write.f:results_ntu_write (l.1-271); util/root_tree/make_root_tree.f + froot.c
// (TTree "h10", nt->Branch(tag, &float, tag) => one Float_t leaf per tag).
#include "simc/output/RootTreeSink.h"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include "simc/core/Registry.h"

namespace simc {

RootTreeSink::RootTreeSink() = default;
RootTreeSink::~RootTreeSink() {
    if (tfile_) tfile_->Close();
    if (worker_ && !file_.empty()) {
        std::error_code ec;
        std::filesystem::remove(file_, ec);   // (a worker that was never merged: error exit)
    }
}

void RootTreeSink::configure(const Json& config, const SimcParameters&) {
    if (config.is_object()) {
        file_ = config.value("file", file_);
        tree_name_ = config.value("tree", tree_name_);
    }
}

void RootTreeSink::init(const StageContext&) {}

void RootTreeSink::beginRun(const RunSetup& setup) {
    setup_ = setup;
    const auto& cols = NtupleColumnsFor(setup);   // throws for the layouts not implemented yet
    if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".root";
    tfile_.reset(TFile::Open(file_.c_str(), "RECREATE"));
    if (!tfile_ || tfile_->IsZombie()) throw std::runtime_error("RootTreeSink: cannot create " + file_);
    tree_ = new TTree(tree_name_.c_str(), tree_name_.c_str());   // froot.c: new TTree("h10","h10")
    ntu_.assign(cols.size(), 0.0);
    out_.assign(cols.size(), 0.0f);
    for (size_t i = 0; i < cols.size(); ++i) tree_->Branch(cols[i].name, &out_[i], cols[i].name);
    entries_ = 0;
}

void RootTreeSink::fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                        bool success) {
    if (!success || !tree_) return;   // results_write.f:34
    FillNtuple(setup_, main, vertex, orig, recon, ntu_);
    for (size_t i = 0; i < ntu_.size(); ++i) out_[i] = static_cast<float>(ntu_[i]);   // make_root_tree real*4
    tree_->Fill();
    ++entries_;
}

void RootTreeSink::endRun(const RunSummary&) { closeFile(); }

void RootTreeSink::closeFile() {
    if (!tfile_) return;
    tfile_->cd();
    tree_->Write();
    tfile_->Close();
    tfile_.reset();
    tree_ = nullptr;
}

// P15.4: parallel event loop (IParallelSink).
std::unique_ptr<IOutputSink> RootTreeSink::makeWorkerSink(int worker, int) const {
    ROOT::EnableThreadSafety();   // every worker fills its own TTree in its own TFile, in its own thread
    auto w = std::make_unique<RootTreeSink>();
    w->file_ = file_ + ".w" + std::to_string(worker) + ".tmp";
    w->tree_name_ = tree_name_;
    w->worker_ = true;
    return w;
}

void RootTreeSink::mergeWorkerSink(IOutputSink& sink) {
    auto& w = dynamic_cast<RootTreeSink&>(sink);
    if (!tree_) throw std::logic_error("RootTreeSink::mergeWorkerSink: the run's tree is not open");
    w.closeFile();
    {
        std::unique_ptr<TFile> f(TFile::Open(w.file_.c_str(), "READ"));
        if (!f || f->IsZombie()) throw std::runtime_error("RootTreeSink: cannot read the worker file " + w.file_);
        TTree* t = f->Get<TTree>(w.tree_name_.c_str());
        if (!t) throw std::runtime_error("RootTreeSink: no tree '" + w.tree_name_ + "' in " + w.file_);
        // the worker's columns are this tree's columns (same layout): read them straight into out_ and Fill
        const auto& cols = NtupleColumnsFor(setup_);
        for (size_t i = 0; i < cols.size(); ++i) t->SetBranchAddress(cols[i].name, &out_[i]);
        const Long64_t n = t->GetEntries();
        for (Long64_t e = 0; e < n; ++e) {
            t->GetEntry(e);
            tree_->Fill();
            ++entries_;
        }
        t->ResetBranchAddresses();
        f->Close();
    }
    std::error_code ec;
    std::filesystem::remove(w.file_, ec);
    w.file_.clear();
}

}  // namespace simc

SIMC_REGISTER(simc::IOutputSink, "root_tree", simc::RootTreeSink)
