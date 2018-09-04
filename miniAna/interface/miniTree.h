#ifndef MINITREE_H
#define MINITREE_H

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <boost/math/special_functions/round.hpp>
#include <memory>
#include <set>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "trauma/miniAna/interface/Event.h"
#include "trauma/miniAna/interface/physicsObjects.h"
#include "trauma/miniAna/interface/eventSelection.h"
#include "trauma/miniAna/interface/configuration.h"


class miniTree {
  public:
    // Default - so root can load based on a name;
    miniTree(configuration &cmaConfig);

    // Default - so we can clean up;
    virtual ~miniTree();

    // Run once at the start of the job;
    virtual void initialize(TFile& outputFile);

    // Run for every event (in every systematic) that needs saving;
    virtual void saveEvent(const std::map<std::string,double> features);

    // Clear stuff;
    virtual void finalize();


  protected:

    TTree * m_ttree;
    TTree * m_metadataTree;
    configuration * m_config;

    /**** Training branches ****/
    // weights for inputs
    float m_weight;
    float m_nominal_weight;

    // Deep learning features
    unsigned int m_target;

    float m_mu_pt;
    float m_mu_eta;
    float m_mu_phi;
    int m_mu_charge;
    float m_tk_pt;
    float m_tk_eta;
    float m_tk_phi;
    float m_tk_sinheta;
    float m_tk_rinv;
    float m_tk_chi2;
    float m_tk_z0;
    float m_tk_d0;
    int m_tk_charge;
    float m_tk_stubPtConsistency;
    float m_tkmu_deltaR2;

    /**** Metadata ****/
    // which sample has which target value
    // many ROOT files will be merged together to do the training!
    std::string m_name;
    unsigned int m_target_value;
    unsigned int m_nEvents;
};

#endif
