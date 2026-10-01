// src/output/RootTreeSink.cpp
// Fortran: results_write.f:results_ntu_write (l.1-271); util/root_tree/make_root_tree.f + froot.c
// (TTree "h10", nt->Branch(tag, &float, tag) => one Float_t leaf per tag).
#include "simc/output/RootTreeSink.h"

#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "simc/core/Registry.h"

namespace simc {

RootTreeSink::RootTreeSink() = default;
RootTreeSink::~RootTreeSink() {
    if (tfile_) tfile_->Close();
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
    if (!UsesEepLayout(setup.flags))
        throw std::runtime_error("RootTreeSink: only the (e,e'p)/phase-space ntuple layout (NtupleInit.f:296-321) "
                                 "is implemented so far");
    if (file_.empty()) file_ = setup.output_dir + "/" + setup.run_name + ".root";
    tfile_.reset(TFile::Open(file_.c_str(), "RECREATE"));
    if (!tfile_ || tfile_->IsZombie()) throw std::runtime_error("RootTreeSink: cannot create " + file_);
    tree_ = new TTree(tree_name_.c_str(), tree_name_.c_str());   // froot.c: new TTree("h10","h10")
    const auto& cols = EepNtupleColumns();
    for (size_t i = 0; i < cols.size(); ++i) tree_->Branch(cols[i].name, &out_[i], cols[i].name);
    entries_ = 0;
}

void RootTreeSink::fill(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                        bool success) {
    if (!success || !tree_) return;   // results_write.f:34
    FillEepNtuple(setup_, main, vertex, orig, recon, ntu_);
    for (size_t i = 0; i < ntu_.size(); ++i) out_[i] = static_cast<float>(ntu_[i]);   // make_root_tree real*4
    tree_->Fill();
    ++entries_;
}

void RootTreeSink::endRun(const RunSummary&) {
    if (!tfile_) return;
    tfile_->cd();
    tree_->Write();
    tfile_->Close();
    tfile_.reset();
    tree_ = nullptr;
}

}  // namespace simc

SIMC_REGISTER(simc::IOutputSink, "root_tree", simc::RootTreeSink)
