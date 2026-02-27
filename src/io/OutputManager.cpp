// src/io/OutputManager.cpp
// Implementation of ROOT output management
// UPDATED: Added intelligent histogram binning for xbj, W, Q2, and physics quantities

#include "simc/io/OutputManager.h"
#include "simc/core/SimcConstants.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <TTree.h>
#include <TDirectory.h>
#include <TGaxis.h>

namespace simc {

// ============================================================================
// Constructor and Destructor
// ============================================================================
OutputManager::OutputManager(const std::string& filename, const std::string& option)
    : filename_(filename) {
    file_ = std::make_unique<TFile>(filename.c_str(), option.c_str());
    
    if (!file_ || file_->IsZombie()) {
        throw std::runtime_error("Cannot open ROOT file: " + filename);
    }
}

OutputManager::~OutputManager() {
    if (initialized_ && file_ && file_->IsOpen()) {
        Finalize();
    }
}

// ============================================================================
// Initialization
// ============================================================================
bool OutputManager::Initialize() {
    if (initialized_) {
        std::cerr << "Warning: OutputManager already initialized" << std::endl;
        return true;
    }
    
    if (!file_ || !file_->IsOpen()) {
        std::cerr << "Error: ROOT file not open" << std::endl;
        return false;
    }
    
    file_->cd();
    
    // Set global axis formatting to avoid excessive precision
    TGaxis::SetMaxDigits(4);
    
    // Create trees
    SetupEventTree();
    SetupGeneratedTree();
    SetupStatsTree();
    
    // Create histograms
    CreateHistograms();
    
    initialized_ = true;
    return true;
}

// ============================================================================
// Tree Setup
// ============================================================================
void OutputManager::SetupEventTree() {
    tree_event_ = new TTree("T", "SIMC Monte Carlo Events");
    tree_event_->SetAutoSave(10000);  // Autosave every 10k events
    
    // Beam and kinematics
    tree_event_->Branch("Ein", &event_.Ein);
    tree_event_->Branch("Q2", &event_.Q2);
    tree_event_->Branch("W", &event_.W);
    tree_event_->Branch("nu", &event_.nu);
    tree_event_->Branch("q", &event_.q);
    tree_event_->Branch("xbj", &event_.xbj);
    tree_event_->Branch("epsilon", &event_.epsilon);
    
    // Missing quantities
    tree_event_->Branch("Em", &event_.Em);
    tree_event_->Branch("Pm", &event_.Pm);
    tree_event_->Branch("Pmiss", &event_.Pmiss);
    tree_event_->Branch("Emiss", &event_.Emiss);
    tree_event_->Branch("Pmx", &event_.Pmx);
    tree_event_->Branch("Pmy", &event_.Pmy);
    tree_event_->Branch("Pmz", &event_.Pmz);
    
    // Recoil
    tree_event_->Branch("Mrec", &event_.Mrec);
    tree_event_->Branch("Trec", &event_.Trec);
    
    // Angles
    tree_event_->Branch("theta_pq", &event_.theta_pq);
    tree_event_->Branch("phi_pq", &event_.phi_pq);
    
    // Electron arm
    tree_event_->Branch("e_delta", &event_.e_delta);
    tree_event_->Branch("e_xptar", &event_.e_xptar);
    tree_event_->Branch("e_yptar", &event_.e_yptar);
    tree_event_->Branch("e_z", &event_.e_z);
    tree_event_->Branch("e_theta", &event_.e_theta);
    tree_event_->Branch("e_phi", &event_.e_phi);
    tree_event_->Branch("e_E", &event_.e_E);
    tree_event_->Branch("e_P", &event_.e_P);
    
    // Hadron arm
    tree_event_->Branch("p_delta", &event_.p_delta);
    tree_event_->Branch("p_xptar", &event_.p_xptar);
    tree_event_->Branch("p_yptar", &event_.p_yptar);
    tree_event_->Branch("p_z", &event_.p_z);
    tree_event_->Branch("p_theta", &event_.p_theta);
    tree_event_->Branch("p_phi", &event_.p_phi);
    tree_event_->Branch("p_E", &event_.p_E);
    tree_event_->Branch("p_P", &event_.p_P);
    
    // Weights
    tree_event_->Branch("weight", &main_.weight);
    tree_event_->Branch("sigcc", &main_.sigcc);
    tree_event_->Branch("jacobian", &main_.jacobian);
    
    // Flags
    tree_event_->Branch("success", &main_.success);
    
    // ========================================================================
    // PHASE 5e: Focal Plane Coordinates
    // ========================================================================
    // Electron focal plane
    tree_event_->Branch("FP_e_x", &main_.FP_electron.x);
    tree_event_->Branch("FP_e_y", &main_.FP_electron.y);
    tree_event_->Branch("FP_e_dx", &main_.FP_electron.dx);
    tree_event_->Branch("FP_e_dy", &main_.FP_electron.dy);
    tree_event_->Branch("FP_e_path", &main_.FP_electron.path);
    
    // Hadron focal plane
    tree_event_->Branch("FP_p_x", &main_.FP_hadron.x);
    tree_event_->Branch("FP_p_y", &main_.FP_hadron.y);
    tree_event_->Branch("FP_p_dx", &main_.FP_hadron.dx);
    tree_event_->Branch("FP_p_dy", &main_.FP_hadron.dy);
    tree_event_->Branch("FP_p_path", &main_.FP_hadron.path);
}

void OutputManager::SetupGeneratedTree() {
    tree_generated_ = new TTree("T_gen", "Generated Events (before cuts)");
    tree_generated_->SetAutoSave(10000);
    
    // Store only generated quantities
    tree_generated_->Branch("Ein", &event_.Ein);
    tree_generated_->Branch("Q2", &event_.Q2);
    tree_generated_->Branch("W", &event_.W);
    tree_generated_->Branch("e_delta", &event_.e_delta);
    tree_generated_->Branch("e_xptar", &event_.e_xptar);
    tree_generated_->Branch("e_yptar", &event_.e_yptar);
    tree_generated_->Branch("p_delta", &event_.p_delta);
    tree_generated_->Branch("p_xptar", &event_.p_xptar);
    tree_generated_->Branch("p_yptar", &event_.p_yptar);
    tree_generated_->Branch("Em", &event_.Em);
    tree_generated_->Branch("Pm", &event_.Pm);
}

void OutputManager::SetupStatsTree() {
    tree_stats_ = new TTree("T_stats", "Run Statistics");
    
    tree_stats_->Branch("nevents_tried", &nevents_tried_);
    tree_stats_->Branch("nevents_generated", &nevents_generated_);
    tree_stats_->Branch("nevents_contribute", &nevents_contribute_);
    tree_stats_->Branch("nevents_written", &nevents_written_);
}

// ============================================================================
// Event Filling
// ============================================================================
bool OutputManager::FillEvent(const SimcEvent& evt, const MainEvent& main) {
    if (!initialized_) {
        std::cerr << "Error: OutputManager not initialized" << std::endl;
        return false;
    }
    
    event_ = evt;
    main_ = main;
    
    tree_event_->Fill();
    ++nevents_written_;
    
    return true;
}

bool OutputManager::FillGenerated(const SimcEvent& evt, double weight) {
    if (!initialized_) {
        return false;
    }
    
    event_ = evt;
    main_.weight = weight;
    
    tree_generated_->Fill();
    ++nevents_generated_;
    
    return true;
}

bool OutputManager::FillReconstructed(const SimcEvent& evt, const MainEvent& main) {
    if (main.success) {
        ++nevents_contribute_;
    }
    
    return FillEvent(evt, main);
}

// ============================================================================
// Histogram Creation - WITH INTELLIGENT BINNING
// ============================================================================
void OutputManager::CreateHistograms() {
    file_->cd();
    
    // Create directory for histograms
    TDirectory* hdir = file_->mkdir("histograms");
    hdir->cd();
    
    // ========================================================================
    // Physics Kinematics - WITH INTELLIGENT BINNING
    // ========================================================================
    
    // xbj: For elastic H(e,e'p), should be ~1.0
    // Use COARSE bins (0.005 width) to avoid floating-point noise artifacts
    hist_gen_1d_["xbj"] = CreateHist1D("h_gen_xbj", 
        "Bjorken x;x_{Bj};Counts", 
        100, 0.0, 2.0);  // Full range
    
    hist_gen_1d_["xbj_zoom"] = CreateHist1D("h_gen_xbj_zoom", 
        "Bjorken x (elastic region);x_{Bj};Counts", 
        20, 0.95, 1.05);  // 0.005 bin width - prevents multiple peaks at "1"
    
    // W: For elastic, should be ~938.272 MeV
    // Use 0.1 MeV bins to show peak structure without over-precision
    hist_gen_1d_["W"] = CreateHist1D("h_gen_W", 
        "Invariant Mass;W (MeV);Counts", 
        200, 800, 2000);  // Full range
    
    hist_gen_1d_["W_elastic"] = CreateHist1D("h_gen_W_elastic", 
        "Invariant Mass (elastic peak);W (MeV);Counts", 
        40, 936, 940);  // 0.1 MeV bins - shows peak without ridiculous precision
    
    // Q²: Wide dynamic range
    hist_gen_1d_["Q2"] = CreateHist1D("h_gen_Q2", 
        "Four-Momentum Transfer;Q^{2} (MeV^{2});Counts", 
        100, 0, 1e7);  // 0 to 10 GeV² in MeV²
    
    hist_gen_1d_["Q2_GeV2"] = CreateHist1D("h_gen_Q2_GeV2", 
        "Four-Momentum Transfer;Q^{2} (GeV^{2});Counts", 
        100, 0, 10);  // In GeV² units
    
    // nu (energy transfer)
    hist_gen_1d_["nu"] = CreateHist1D("h_gen_nu", 
        "Energy Transfer;#nu (MeV);Counts", 
        100, 0, 5000);
    
    // epsilon (virtual photon polarization)
    hist_gen_1d_["epsilon"] = CreateHist1D("h_gen_epsilon", 
        "Virtual Photon Polarization;#epsilon;Counts", 
        100, 0, 1);
    
    // ========================================================================
    // Energy Conservation Check
    // ========================================================================
    
    // E_in - E_e' - E_p + M_p should be ~0 for elastic
    hist_gen_1d_["Econs"] = CreateHist1D("h_gen_Econs", 
        "Energy Conservation;E_{in} - E_{e'} - E_{p} + M_{p} (MeV);Counts", 
        100, -10, 10);
    
    // ========================================================================
    // Generated - Electron arm
    // ========================================================================
    
    hist_gen_1d_["e_delta"] = CreateHist1D("h_gen_e_delta", 
        "Generated Electron #delta;#delta (%);Counts", 
        100, -15, 15);
    hist_gen_1d_["e_xptar"] = CreateHist1D("h_gen_e_xptar", 
        "Generated Electron x'_{tar};x'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_gen_1d_["e_yptar"] = CreateHist1D("h_gen_e_yptar", 
        "Generated Electron y'_{tar};y'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_gen_1d_["e_theta"] = CreateHist1D("h_gen_e_theta", 
        "Generated Electron Angle;#theta_{e} (rad);Counts", 
        100, 0, 0.5);
    hist_gen_1d_["e_E"] = CreateHist1D("h_gen_e_E", 
        "Generated Electron Energy;E_{e'} (MeV);Counts", 
        100, 0, 12000);
    
    // ========================================================================
    // Generated - Hadron arm
    // ========================================================================
    
    hist_gen_1d_["p_delta"] = CreateHist1D("h_gen_p_delta", 
        "Generated Hadron #delta;#delta (%);Counts", 
        100, -20, 20);
    hist_gen_1d_["p_xptar"] = CreateHist1D("h_gen_p_xptar", 
        "Generated Hadron x'_{tar};x'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_gen_1d_["p_yptar"] = CreateHist1D("h_gen_p_yptar", 
        "Generated Hadron y'_{tar};y'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_gen_1d_["p_theta"] = CreateHist1D("h_gen_p_theta", 
        "Generated Proton Angle;#theta_{p} (rad);Counts", 
        100, 0, 1.5);
    hist_gen_1d_["p_E"] = CreateHist1D("h_gen_p_E", 
        "Generated Proton Energy;E_{p} (MeV);Counts", 
        100, 0, 8000);
    
    // ========================================================================
    // Missing Momentum
    // ========================================================================
    
    hist_gen_1d_["Pmiss"] = CreateHist1D("h_gen_Pmiss", 
        "Missing Momentum;P_{miss} (MeV/c);Counts", 
        100, 0, 500);
    hist_gen_1d_["Em"] = CreateHist1D("h_gen_Em", 
        "Missing Energy;E_{miss} (MeV);Counts", 
        100, -100, 200);
    
    // ========================================================================
    // Reconstructed - Electron arm
    // ========================================================================
    
    hist_rec_1d_["e_delta"] = CreateHist1D("h_rec_e_delta", 
        "Reconstructed Electron #delta;#delta (%);Counts", 
        100, -15, 15);
    hist_rec_1d_["e_xptar"] = CreateHist1D("h_rec_e_xptar", 
        "Reconstructed Electron x'_{tar};x'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_rec_1d_["e_yptar"] = CreateHist1D("h_rec_e_yptar", 
        "Reconstructed Electron y'_{tar};y'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_rec_1d_["e_theta"] = CreateHist1D("h_rec_e_theta", 
        "Reconstructed Electron Angle;#theta_{e} (rad);Counts", 
        100, 0, 0.5);
    hist_rec_1d_["e_E"] = CreateHist1D("h_rec_e_E", 
        "Reconstructed Electron Energy;E_{e'} (MeV);Counts", 
        100, 0, 12000);
    
    // ========================================================================
    // Reconstructed - Hadron arm
    // ========================================================================
    
    hist_rec_1d_["p_delta"] = CreateHist1D("h_rec_p_delta", 
        "Reconstructed Hadron #delta;#delta (%);Counts", 
        100, -20, 20);
    hist_rec_1d_["p_xptar"] = CreateHist1D("h_rec_p_xptar", 
        "Reconstructed Hadron x'_{tar};x'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_rec_1d_["p_yptar"] = CreateHist1D("h_rec_p_yptar", 
        "Reconstructed Hadron y'_{tar};y'_{tar} (rad);Counts", 
        100, -0.1, 0.1);
    hist_rec_1d_["p_theta"] = CreateHist1D("h_rec_p_theta", 
        "Reconstructed Proton Angle;#theta_{p} (rad);Counts", 
        100, 0, 1.5);
    hist_rec_1d_["p_E"] = CreateHist1D("h_rec_p_E", 
        "Reconstructed Proton Energy;E_{p} (MeV);Counts", 
        100, 0, 8000);
    
    // ========================================================================
    // Reconstructed - Missing momentum
    // ========================================================================
    
    hist_rec_1d_["Pmiss"] = CreateHist1D("h_rec_Pmiss", 
        "Reconstructed Missing Momentum;P_{miss} (MeV/c);Counts", 
        100, 0, 500);
    hist_rec_1d_["Em"] = CreateHist1D("h_rec_Em", 
        "Reconstructed Missing Energy;E_{miss} (MeV);Counts", 
        100, -100, 200);
    
    // ========================================================================
    // 2D histograms - Generated
    // ========================================================================
    
    hist_gen_2d_["Em_vs_Pm"] = CreateHist2D("h_gen_Em_vs_Pm", 
        "E_{miss} vs P_{miss};P_{miss} (MeV/c);E_{miss} (MeV)",
        100, 0, 500, 100, -100, 200);
    
    hist_gen_2d_["Q2_vs_W"] = CreateHist2D("h_gen_Q2_vs_W", 
        "Q^{2} vs W;W (MeV);Q^{2} (GeV^{2})",
        100, 800, 2000, 100, 0, 5);
    
    hist_gen_2d_["xbj_vs_Q2"] = CreateHist2D("h_gen_xbj_vs_Q2", 
        "x_{Bj} vs Q^{2};Q^{2} (GeV^{2});x_{Bj}",
        100, 0, 10, 100, 0, 2);
    
    // ========================================================================
    // 2D histograms - Reconstructed
    // ========================================================================
    
    hist_rec_2d_["Em_vs_Pm"] = CreateHist2D("h_rec_Em_vs_Pm", 
        "E_{miss} vs P_{miss} (Reconstructed);P_{miss} (MeV/c);E_{miss} (MeV)",
        100, 0, 500, 100, -100, 200);
    
    hist_rec_2d_["Q2_vs_W"] = CreateHist2D("h_rec_Q2_vs_W", 
        "Q^{2} vs W (Reconstructed);W (MeV);Q^{2} (GeV^{2})",
        100, 800, 2000, 100, 0, 5);
    
    file_->cd();
}

TH1D* OutputManager::CreateHist1D(const std::string& name, const std::string& title,
                                 int nbins, double xmin, double xmax) {
    return new TH1D(name.c_str(), title.c_str(), nbins, xmin, xmax);
}

TH2D* OutputManager::CreateHist2D(const std::string& name, const std::string& title,
                                 int nbinsx, double xmin, double xmax,
                                 int nbinsy, double ymin, double ymax) {
    return new TH2D(name.c_str(), title.c_str(), 
                    nbinsx, xmin, xmax, 
                    nbinsy, ymin, ymax);
}

// ============================================================================
// Histogram Filling - UPDATED TO INCLUDE NEW HISTOGRAMS
// ============================================================================
void OutputManager::FillHistograms(const SimcEvent& evt, double weight, bool generated) {
    auto& hist_map = generated ? hist_gen_1d_ : hist_rec_1d_;
    auto& hist_map_2d = generated ? hist_gen_2d_ : hist_rec_2d_;
    
    // Fill physics kinematics (generated only)
    if (generated) {
        if (hist_map.count("xbj")) hist_map["xbj"]->Fill(evt.xbj, weight);
        if (hist_map.count("xbj_zoom")) hist_map["xbj_zoom"]->Fill(evt.xbj, weight);
        if (hist_map.count("W")) hist_map["W"]->Fill(evt.W, weight);
        if (hist_map.count("W_elastic")) hist_map["W_elastic"]->Fill(evt.W, weight);
        if (hist_map.count("Q2")) hist_map["Q2"]->Fill(evt.Q2, weight);
        if (hist_map.count("Q2_GeV2")) hist_map["Q2_GeV2"]->Fill(evt.Q2/1e6, weight);
        if (hist_map.count("nu")) hist_map["nu"]->Fill(evt.nu, weight);
        if (hist_map.count("epsilon")) hist_map["epsilon"]->Fill(evt.epsilon, weight);
        
        // Energy conservation: E_in - E_e' - E_p + M_p
        double Econs = evt.Ein - evt.e_E - evt.p_E + constants::Mp;  // Mp from SimcConstants
        if (hist_map.count("Econs")) hist_map["Econs"]->Fill(Econs, weight);
        
        // Angles and energies
        if (hist_map.count("e_theta")) hist_map["e_theta"]->Fill(evt.e_theta, weight);
        if (hist_map.count("e_E")) hist_map["e_E"]->Fill(evt.e_E, weight);
        if (hist_map.count("p_theta")) hist_map["p_theta"]->Fill(evt.p_theta, weight);
        if (hist_map.count("p_E")) hist_map["p_E"]->Fill(evt.p_E, weight);
        
        // Missing momentum
        if (hist_map.count("Pmiss")) hist_map["Pmiss"]->Fill(evt.Pmiss, weight);
        if (hist_map.count("Em")) hist_map["Em"]->Fill(evt.Em, weight);
        
        // 2D histograms (generated only)
        if (hist_map_2d.count("xbj_vs_Q2")) {
            hist_map_2d["xbj_vs_Q2"]->Fill(evt.Q2/1e6, evt.xbj, weight);
        }
    }
    
    // Fill spectrometer quantities (both generated and reconstructed)
    if (hist_map.count("e_delta")) hist_map["e_delta"]->Fill(evt.e_delta, weight);
    if (hist_map.count("e_xptar")) hist_map["e_xptar"]->Fill(evt.e_xptar, weight);
    if (hist_map.count("e_yptar")) hist_map["e_yptar"]->Fill(evt.e_yptar, weight);
    if (hist_map.count("p_delta")) hist_map["p_delta"]->Fill(evt.p_delta, weight);
    if (hist_map.count("p_xptar")) hist_map["p_xptar"]->Fill(evt.p_xptar, weight);
    if (hist_map.count("p_yptar")) hist_map["p_yptar"]->Fill(evt.p_yptar, weight);
    
    // Fill angles and energies (both generated and reconstructed)
    if (hist_map.count("e_theta")) hist_map["e_theta"]->Fill(evt.e_theta, weight);
    if (hist_map.count("e_E")) hist_map["e_E"]->Fill(evt.e_E, weight);
    if (hist_map.count("p_theta")) hist_map["p_theta"]->Fill(evt.p_theta, weight);
    if (hist_map.count("p_E")) hist_map["p_E"]->Fill(evt.p_E, weight);
    
    // Fill missing momentum (both generated and reconstructed)
    if (hist_map.count("Pmiss")) hist_map["Pmiss"]->Fill(evt.Pmiss, weight);
    if (hist_map.count("Em")) hist_map["Em"]->Fill(evt.Em, weight);
    
    // Fill 2D histograms (both generated and reconstructed)
    if (hist_map_2d.count("Em_vs_Pm")) {
        hist_map_2d["Em_vs_Pm"]->Fill(evt.Pm, evt.Em, weight);
    }
    if (hist_map_2d.count("Q2_vs_W")) {
        hist_map_2d["Q2_vs_W"]->Fill(evt.W, evt.Q2/1e6, weight);
    }
}

void OutputManager::WriteHistograms() {
    file_->cd("histograms");
    
    for (auto& [name, hist] : hist_gen_1d_) {
        hist->Write();
    }
    for (auto& [name, hist] : hist_gen_2d_) {
        hist->Write();
    }
    for (auto& [name, hist] : hist_rec_1d_) {
        hist->Write();
    }
    for (auto& [name, hist] : hist_rec_2d_) {
        hist->Write();
    }
    
    file_->cd();
}

// ============================================================================
// Statistics
// ============================================================================
void OutputManager::SetStatistics(long long ntried, long long ngenerated, long long ncontribute) {
    nevents_tried_ = ntried;
    nevents_generated_ = ngenerated;
    nevents_contribute_ = ncontribute;
}

// ============================================================================
// Finalization
// ============================================================================
bool OutputManager::Finalize() {
    if (!initialized_) {
        return false;
    }
    
    if (!file_ || !file_->IsOpen()) {
        std::cerr << "Error: ROOT file not open" << std::endl;
        return false;
    }
    
    file_->cd();
    
    // Fill statistics tree
    tree_stats_->Fill();
    
    // Write all trees
    if (tree_event_) tree_event_->Write();
    if (tree_generated_) tree_generated_->Write();
    if (tree_stats_) tree_stats_->Write();
    
    // Write histograms
    WriteHistograms();
    
    // Print summary
    std::cout << "\n========== SIMC Output Summary ==========" << std::endl;
    std::cout << "Events tried:       " << nevents_tried_ << std::endl;
    std::cout << "Events generated:   " << nevents_generated_ << std::endl;
    std::cout << "Events contribute:  " << nevents_contribute_ << std::endl;
    std::cout << "Events written:     " << nevents_written_ << std::endl;
    std::cout << "Output file:        " << filename_ << std::endl;
    std::cout << "=========================================\n" << std::endl;
    
    // Close file
    file_->Close();
    
    initialized_ = false;
    return true;
}

} // namespace simc
