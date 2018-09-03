#ifndef EVENT_H
#define EVENT_H

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "trauma/miniAna/interface/physicsObjects.h"
#include "trauma/miniAna/interface/configuration.h"
#include "trauma/miniAna/interface/deepLearning.h"


// Event Class
class Event {
  public:
    // Constructor
    Event( TTreeReader &myReader, configuration &cmaConfig);
    Event( const Event &obj);

    // Destructor
    virtual ~Event();

    // Execute the event (load information and setup objects)
    virtual void execute(Long64_t entry);
    virtual void updateEntry(Long64_t entry);
    bool isValidRecoEntry() const {return (m_entry > (long long)-1);}

    // Setup physics information
    void initialize_tracks();
    void initialize_muons();
    void initialize_weights();

    // Clear stuff;
    virtual void finalize();
    virtual void clear();

    // Get physics object(s) information
    std::vector<Muon> bmtf_muons() const {return m_bmtf_muons;}
    std::vector<Muon> omtf_muons() const {return m_omtf_muons;}
    std::vector<Muon> emtf_muons() const {return m_emtf_muons;}
    std::vector<Track> tracks() const {return m_tracks;}
    std::vector<TkMu> trackmuons() const {return m_tkmuons;}

    // metadata
    long long entry() const {return m_entry;}
    virtual std::string treeName() const {return m_treeName;}

    // functions to calculate things
    void deepLearningPrediction();

    // Get weights
    virtual float nominal_weight() const {return m_nominal_weight;}

  protected:

    // general information
    configuration *m_config;
    TTreeReader &m_ttree;
    TTreeReader m_truth_tree;
    std::string m_treeName;
    bool m_isMC;
    long long m_entry;
    long long m_truth_entry;

    // event weight information
    double m_nominal_weight;
    double m_xsection;
    double m_kfactor;
    double m_sumOfWeights;
    double m_LUMI;

    // physics object information
    std::vector<Muon> m_bmtf_muons;
    std::vector<Muon> m_omtf_muons;
    std::vector<Muon> m_emtf_muons;
    std::vector<Track> m_tracks;
    std::vector<TkMu> m_tkmuons;

    deepLearning* m_deepLearningTool;            // tool to perform deep learning
    bool m_DNNtraining;
    bool m_DNNinference;

    // TTree variables
    // *************
    // Muon info
    TTreeReaderValue<std::vector<float>> * m_mu_bmtf_pt;
    TTreeReaderValue<std::vector<float>> * m_mu_bmtf_eta;
    TTreeReaderValue<std::vector<float>> * m_mu_bmtf_phi;
    TTreeReaderValue<std::vector<float>> * m_mu_bmtf_e;
    TTreeReaderValue<std::vector<float>> * m_mu_bmtf_charge;

    TTreeReaderValue<std::vector<float>> * m_mu_omtf_pt;
    TTreeReaderValue<std::vector<float>> * m_mu_omtf_eta;
    TTreeReaderValue<std::vector<float>> * m_mu_omtf_phi;
    TTreeReaderValue<std::vector<float>> * m_mu_omtf_e;
    TTreeReaderValue<std::vector<float>> * m_mu_omtf_charge;

    TTreeReaderValue<std::vector<float>> * m_mu_emtf_pt;
    TTreeReaderValue<std::vector<float>> * m_mu_emtf_eta;
    TTreeReaderValue<std::vector<float>> * m_mu_emtf_phi;
    TTreeReaderValue<std::vector<float>> * m_mu_emtf_e;
    TTreeReaderValue<std::vector<float>> * m_mu_emtf_charge;

    // TTTrack info
    TTreeReaderValue<std::vector<float>> * m_track_pt;
    TTreeReaderValue<std::vector<float>> * m_track_eta;
    TTreeReaderValue<std::vector<float>> * m_track_phi;
    TTreeReaderValue<std::vector<float>> * m_track_d0;
    TTreeReaderValue<std::vector<float>> * m_track_z0;
    TTreeReaderValue<std::vector<float>> * m_track_chi2;
    TTreeReaderValue<std::vector<int>> * m_track_charge;
    TTreeReaderValue<std::vector<int>> * m_track_nstub;
    TTreeReaderValue<std::vector<int>> * m_track_genuine;
    TTreeReaderValue<std::vector<int>> * m_track_loose;
    TTreeReaderValue<std::vector<int>> * m_track_unknown;
    TTreeReaderValue<std::vector<int>> * m_track_combinatoric;
    TTreeReaderValue<std::vector<int>> * m_track_fake;
    TTreeReaderValue<std::vector<int>> * m_track_matchtp_pdgid;
    TTreeReaderValue<std::vector<float>> * m_track_matchtp_pt;
    TTreeReaderValue<std::vector<float>> * m_track_matchtp_eta;
    TTreeReaderValue<std::vector<float>> * m_track_matchtp_phi;
    TTreeReaderValue<std::vector<float>> * m_track_matchtp_z0;
    TTreeReaderValue<std::vector<float>> * m_track_matchtp_dxy;
};

#endif
