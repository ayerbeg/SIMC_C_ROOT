// include/simc/OutputManager.h
// ROOT-based output management for SIMC Monte Carlo

#ifndef SIMC_OUTPUT_MANAGER_H
#define SIMC_OUTPUT_MANAGER_H

#include "simc/core/SimcEvent.h"
#include <string>
#include <memory>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>

namespace simc {

/**
 * @class OutputManager
 * @brief Manages ROOT file output including TTrees and histograms
 * 
 * This class handles all output operations for SIMC, including:
 * - Creating and filling ROOT TTrees with event data
 * - Managing diagnostic histograms
 * - Writing summary information
 * - Proper file closure and cleanup
 * 
 * The TTree structure is optimized for analysis with flat branches
 * rather than complex nested objects.
 * 
 * Example usage:
 * @code
 * OutputManager output("simc_output.root");
 * output.Initialize();
 * 
 * for (int i = 0; i < nevents; ++i) {
 *     SimcEvent evt;
 *     MainEvent main;
 *     // ... generate event ...
 *     output.FillEvent(evt, main);
 * }
 * 
 * output.Finalize();
 * @endcode
 */
class OutputManager {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Constructor with output filename
     * @param filename Output ROOT filename
     * @param option File open option ("RECREATE" or "UPDATE")
     */
    explicit OutputManager(const std::string& filename, 
                          const std::string& option = "RECREATE");
    
    /**
     * @brief Destructor - ensures proper file closure
     */
    ~OutputManager();
    
    // ========================================================================
    // Initialization and Finalization
    // ========================================================================
    
    /**
     * @brief Initialize output system (create trees and histograms)
     * @return true if successful
     */
    bool Initialize();
    
    /**
     * @brief Finalize output (write and close file)
     * @return true if successful
     */
    bool Finalize();
    
    /**
     * @brief Check if initialized
     */
    bool IsInitialized() const { return initialized_; }
    
    // ========================================================================
    // Event Output
    // ========================================================================
    
    /**
     * @brief Fill event tree with generated event
     * @param evt Event data
     * @param main Additional event information
     * @return true if successful
     */
    bool FillEvent(const SimcEvent& evt, const MainEvent& main);
    
    /**
     * @brief Fill generated (initial) event tree
     * @param evt Event data
     * @return true if successful
     */
    bool FillGenerated(const SimcEvent& evt, double weight = 1.0);
    
    /**
     * @brief Fill reconstructed event tree
     * @param evt Event data
     * @param main Additional event information
     * @return true if successful
     */
    bool FillReconstructed(const SimcEvent& evt, const MainEvent& main);
    
    // ========================================================================
    // Histogram Management
    // ========================================================================
    
    /**
     * @brief Create standard diagnostic histograms
     */
    void CreateHistograms();
    
    /**
     * @brief Fill diagnostic histograms
     * @param evt Event data
     * @param weight Event weight
     * @param generated True if generated quantities, false if reconstructed
     */
    void FillHistograms(const SimcEvent& evt, double weight, bool generated);
    
    /**
     * @brief Write all histograms to file
     */
    void WriteHistograms();
    
    // ========================================================================
    // Statistics and Counters
    // ========================================================================
    
    /**
     * @brief Get number of events written
     */
    long long GetNEventsWritten() const { return nevents_written_; }
    
    /**
     * @brief Get number of events tried
     */
    long long GetNEventsTried() const { return nevents_tried_; }
    
    /**
     * @brief Increment tried counter
     */
    void IncrementTried() { ++nevents_tried_; }
    
    /**
     * @brief Set simulation statistics
     */
    void SetStatistics(long long ntried, long long ngenerated, long long ncontribute);
    
    // ========================================================================
    // Summary Information
    // ========================================================================
    
    /**
     * @brief Write summary information to file
     * @param info Map of parameter name to value
     */
    void WriteSummary(const std::map<std::string, double>& info);
    
    /**
     * @brief Write configuration to file
     * @param config_string Configuration as string
     */
    void WriteConfiguration(const std::string& config_string);
    
    // ========================================================================
    // File Access
    // ========================================================================
    
    /**
     * @brief Get ROOT file pointer
     */
    TFile* GetFile() { return file_.get(); }
    
    /**
     * @brief Get main event tree
     */
    TTree* GetEventTree() { return tree_event_; }
    
private:
    // ========================================================================
    // Private Members
    // ========================================================================
    
    // File management
    std::unique_ptr<TFile> file_;           ///< ROOT output file
    std::string filename_;                  ///< Output filename
    bool initialized_{false};               ///< Initialization flag
    
    // Trees
    TTree* tree_event_{nullptr};            ///< Main event tree
    TTree* tree_generated_{nullptr};        ///< Generated events only
    TTree* tree_stats_{nullptr};            ///< Statistics tree
    
    // Event data (for TTree branches)
    SimcEvent event_;                       ///< Current event
    MainEvent main_;                        ///< Additional event info
    
    // Counters
    long long nevents_written_{0};          ///< Events written to tree
    long long nevents_tried_{0};            ///< Events attempted
    long long nevents_generated_{0};        ///< Events generated
    long long nevents_contribute_{0};       ///< Events contributing
    
    // Histograms - Generated quantities
    std::map<std::string, TH1D*> hist_gen_1d_;
    std::map<std::string, TH2D*> hist_gen_2d_;
    
    // Histograms - Reconstructed quantities
    std::map<std::string, TH1D*> hist_rec_1d_;
    std::map<std::string, TH2D*> hist_rec_2d_;
    
    // ========================================================================
    // Private Methods
    // ========================================================================
    
    /**
     * @brief Setup event tree branches
     */
    void SetupEventTree();
    
    /**
     * @brief Setup generated event tree branches
     */
    void SetupGeneratedTree();
    
    /**
     * @brief Setup statistics tree
     */
    void SetupStatsTree();
    
    /**
     * @brief Create 1D histogram
     */
    TH1D* CreateHist1D(const std::string& name, const std::string& title,
                       int nbins, double xmin, double xmax);
    
    /**
     * @brief Create 2D histogram
     */
    TH2D* CreateHist2D(const std::string& name, const std::string& title,
                       int nbinsx, double xmin, double xmax,
                       int nbinsy, double ymin, double ymax);
    
    /**
     * @brief Get or create 1D histogram
     */
    TH1D* GetHist1D(const std::string& name, bool generated);
    
    /**
     * @brief Get or create 2D histogram
     */
    TH2D* GetHist2D(const std::string& name, bool generated);
};

/**
 * @class HistogramSet
 * @brief Convenient container for related histograms
 * 
 * Groups related histograms (e.g., all electron arm quantities)
 * for easier management.
 */
class HistogramSet {
public:
    /**
     * @brief Constructor
     * @param name Set name (prefix for histograms)
     */
    explicit HistogramSet(const std::string& name);
    
    /**
     * @brief Add 1D histogram
     */
    void Add1D(const std::string& var, const std::string& title,
               int nbins, double xmin, double xmax);
    
    /**
     * @brief Add 2D histogram
     */
    void Add2D(const std::string& varx, const std::string& vary,
               const std::string& title,
               int nbinsx, double xmin, double xmax,
               int nbinsy, double ymin, double ymax);
    
    /**
     * @brief Fill 1D histogram
     */
    void Fill1D(const std::string& var, double value, double weight = 1.0);
    
    /**
     * @brief Fill 2D histogram
     */
    void Fill2D(const std::string& var, double x, double y, double weight = 1.0);
    
    /**
     * @brief Write all histograms to current directory
     */
    void Write();
    
private:
    std::string name_;
    std::map<std::string, TH1D*> hist_1d_;
    std::map<std::string, TH2D*> hist_2d_;
};

} // namespace simc

#endif // SIMC_OUTPUT_MANAGER_H
