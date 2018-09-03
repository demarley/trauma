/*
Created:         3 September 2018
Last Updated:    3 September 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Create and fill TTree for ML
*/
#include "trauma/miniAna/interface/miniTree.h"


miniTree::miniTree(configuration &cmaConfig) : 
  m_config(&cmaConfig){}

miniTree::~miniTree() {}


void miniTree::initialize(TFile& outputFile) {
    /*
       Setup the new tree 
       Contains features for the NN
       --  No vector<T> stored in tree: completely flat!
    */
    outputFile.cd();                                     // move to output file
    m_ttree        = new TTree("features", "features");  // Tree contains features for the NN
    m_metadataTree = new TTree("metadata","metadata");   // Tree contains metadata

    /**** Setup new branches here ****/
    // Weights
    m_ttree->Branch( "weight",   &m_weight,   "weight/F" );
    m_ttree->Branch( "nominal_weight", &m_nominal_weight, "nominal_weight/F" );

    // Features
    m_ttree->Branch( "target", &m_target, "target/I" );  // target value (.e.g, 0 or 1)

    m_ttree->Branch( "mu_pt",     &m_mu_pt,  "mu_pt/F" );
    m_ttree->Branch( "mu_eta",    &m_mu_eta, "mu_eta/F" );
    m_ttree->Branch( "mu_phi",    &m_mu_phi, "mu_phi/F" );
    m_ttree->Branch( "mu_charge", &m_mu_charge, "mu_charge/F" );

    m_ttree->Branch( "tk_pt",      &m_tk_pt,      "tk_pt/F" );
    m_ttree->Branch( "tk_eta",     &m_tk_eta,     "tk_eta/F" );
    m_ttree->Branch( "tk_phi",     &m_tk_phi,     "tk_phi/F" );
    m_ttree->Branch( "tk_sinheta", &m_tk_sinheta, "tk_sinheta/F" );
    m_ttree->Branch( "tk_rinv",    &m_tk_rinv,    "tk_rinv/F" );
    m_ttree->Branch( "tk_chi2",    &m_tk_chi2,    "tk_chi2/F" );
    m_ttree->Branch( "tk_z0",      &m_tk_z0,      "tk_z0/F" );
    m_ttree->Branch( "tk_d0",      &m_tk_d0,      "tk_d0/F" );
    m_ttree->Branch( "tk_charge",  &m_tk_charge,  "tk_charge/F" );
//    m_ttree->Branch( "tk_stubPtConsistency",     &m_tk_stubPtConsistency,  "tk_stubPtConsistency/F" );

    m_ttree->Branch( "tkmu_deltaR2", &m_tkmu_deltaR2,  "tkmu_deltaR2/F" );

    /**** Metadata ****/
    m_metadataTree->Branch( "nEvents", &m_nEvents, "nEvents/I" );

    return;
} // end initialize



void miniTree::saveEvent(const std::map<std::string,double> features) {
    /* Save the ML features to the ttree! */
    cma::DEBUG("MINITREE : Save event ");

    m_weight   = features.at("weight");
    m_nominal_weight = features.at("nominal_weight");

    m_target = features.at("target");

    m_mu_pt  = features.at("mu_pt");
    m_mu_eta = features.at("mu_eta");
    m_mu_phi = features.at("mu_phi");
    m_mu_charge = features.at("mu_charge");

    m_tk_pt      = features.at("tk_pt");
    m_tk_eta     = features.at("tk_eta");
    m_tk_phi     = features.at("tk_phi");
    m_tk_sinheta = features.at("tk_sinheta");
    m_tk_rinv    = features.at("tk_rinv");
    m_tk_chi2    = features.at("tk_chi2");
    m_tk_z0      = features.at("tk_z0");
    m_tk_d0      = features.at("tk_d0");
    m_tk_charge  = features.at("tk_charge");
    //m_tk_stubPtConsistency;

    m_tkmu_deltaR2 = features.at("tkmu_deltaR2");

    /**** Fill the tree ****/
    cma::DEBUG("MINITREE : Fill the tree");
    m_ttree->Fill();

    return;
}


void miniTree::finalize(){
    /* Finalize the class -- fill in the metadata (only need to do this once!) */
    m_nEvents = m_config->NTotalEvents();

    cma::DEBUG("MINITREE : Fill the metadata tree");
    m_metadataTree->Fill();
}

// THE END
