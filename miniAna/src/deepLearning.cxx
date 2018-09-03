/*
Created:         3 September 2018
Last Updated:    3 September 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Tool for performing deep learning tasks
*/
#include "trauma/miniAna/interface/deepLearning.h"


deepLearning::deepLearning( configuration& cmaConfig ) :
  m_config(&cmaConfig),
  m_lwnn(nullptr),
  m_dnnKey(""){
    if (m_config->DNNinference()){
        // Setup lwtnn
        std::ifstream input_cfg = cma::open_file( m_config->dnnFile() );
        lwt::JSONConfig cfg     = lwt::parse_json( input_cfg );
        m_lwnn   = new lwt::LightweightNeuralNetwork(cfg.inputs, cfg.layers, cfg.outputs);
        m_dnnKey = m_config->dnnKey();
    }
  }

deepLearning::~deepLearning() {
    delete m_lwnn;
}


void deepLearning::training(const Track& tk, const Muon& mu){
    /* Prepare inputs for training */
    m_track = tk;
    m_muon  = mu;
    loadFeatures();
    return;
}

void deepLearning::inference(const Track& tk, const Muon& mu){
    /* Obtain results from LWTNN */
    m_track = tk;
    m_muon  = mu;

    loadFeatures();
    m_discriminant = m_lwnn->compute(m_features);

    return;
}


void deepLearning::loadFeatures(){
    /* Calculate DNN features */
    m_features.clear();

    // feature calculations
    // if matched at gen-level, target = 1 else 0
    unsigned int target(0);
    if (fabs(m_track.matchtp_pdgid)==13){
        target = 1;
    }

    m_features["target"] = target;

    // pt,eta,phi
    // sinh(eta),Rinv,chi2,z0,d0,nstubs  // stub pt consistency (?)
    // deltaR2 (DeltaR ** 2)
    m_features["mu_pt"]  = m_muon.p4.Pt();
    m_features["mu_eta"] = m_muon.p4.Eta();
    m_features["mu_phi"] = m_muon.p4.Phi();
    m_features["mu_charge"] = m_muon.charge;

    m_features["tk_sinheta"] = sinh(m_track.p4.Eta());
    m_features["tk_Rinv"]    = 0.0114/m_track.p4.Pt();      // Rinv = (0.3*3.8/100.0)/pT
    m_features["tk_chi2"]    = m_track.chi2;
    m_features["tk_z0"]      = m_track.z0;
    m_features["tk_d0"]      = m_track.d0;
    m_features["tk_charge"]  = m_track.charge;
//    m_features["tk_nstubs"]  = m_track.nstubs;            // accessible in firmware?
//    m_features["tk_stubPtConsistency"] = m_track.stubPtConsistency;  // not saved in first set of samples >:|

    m_features["tkmu_deltaR2"] = pow((m_muon.p4.DeltaR( m_track.p4 )),2);  // accessible in firmware

    m_features["weight"] = 1.;

    cma::DEBUG("EVENT : Set DNN input values ");

    return;
}

double deepLearning::prediction(){
    /* Return the score for the default key */
    return m_discriminant.at(m_dnnKey);
}

double deepLearning::prediction(const std::string& key){
    /* Just return the prediction (after execute!) */
    return m_discriminant.at(key);
}

// THE END //
