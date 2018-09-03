/*
Created:         3 September 2018
Last Updated:    3 September 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Event class
 Contains all the objects (& structs) with event information
*/
#include "trauma/miniAna/interface/Event.h"

// constructor
Event::Event( TTreeReader &myReader, configuration &cmaConfig ) :
  m_config(&cmaConfig),
  m_ttree(myReader),
  m_treeName("SetMe"){
    m_isMC     = m_config->isMC();
    m_treeName = m_ttree.GetTree()->GetName();              // for systematics
    m_DNNtraining  = m_config->DNNtraining();
    m_DNNinference = m_config->DNNinference();

    // ** LOAD BRANCHES FROM TTREE ** //
    m_mu_bmtf_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_bmtf_pt");
    m_mu_bmtf_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_bmtf_eta");
    m_mu_bmtf_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_bmtf_phi");
    m_mu_bmtf_charge = new TTreeReaderValue<std::vector<int>>(m_ttree,"mu_bmtf_q");

    m_mu_omtf_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_omtf_pt");
    m_mu_omtf_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_omtf_eta");
    m_mu_omtf_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_omtf_phi");
    m_mu_omtf_charge = new TTreeReaderValue<std::vector<int>>(m_ttree,"mu_omtf_q");

    m_mu_emtf_pt  = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_emtf_pt");
    m_mu_emtf_eta = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_emtf_eta");
    m_mu_emtf_phi = new TTreeReaderValue<std::vector<float>>(m_ttree,"mu_emtf_phi");
    m_mu_emtf_charge = new TTreeReaderValue<std::vector<int>>(m_ttree,"mu_emtf_q");

    m_track_pt   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_pt");
    m_track_eta  = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_eta");
    m_track_phi  = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_phi");
    m_track_d0   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_z0");
    m_track_z0   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_z0");
    m_track_chi2 = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_chi2");
    m_track_charge  = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_charge");
    m_track_nstub   = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_nstub");
    m_track_genuine = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_genuine");
    m_track_loose   = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_loose");
    m_track_unknown = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_unknown");
    m_track_fake    = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_fake");
    m_track_combinatoric  = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_combinatoric");
    m_track_matchtp_pdgid = new TTreeReaderValue<std::vector<int>>(m_ttree,"trk_matchtp_pdgid");
    m_track_matchtp_pt    = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_matchtp_pt");
    m_track_matchtp_eta   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_matchtp_eta");
    m_track_matchtp_phi   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_matchtp_phi");
    m_track_matchtp_z0    = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_matchtp_z0");
    m_track_matchtp_dxy   = new TTreeReaderValue<std::vector<float>>(m_ttree,"trk_matchtp_dxy");
} // end constructor

Event::~Event() {}

void Event::clear(){
    /* Clear many of the vectors/maps for each event */
    m_bmtf_muons.clear();
    m_omtf_muons.clear();
    m_emtf_muons.clear();

    m_tracks.clear();

    m_tkmuons.clear();

    return;
}


void Event::updateEntry(Long64_t entry){
    /* Update the entry -> update all TTree variables */
    cma::DEBUG("EVENT : Update Entry "+std::to_string(entry) );

    m_entry = entry;

    // make sure the entry exists/is valid
    if (isValidRecoEntry())
        m_ttree.SetEntry(m_entry);
    else
        cma::ERROR("EVENT : Invalid Reco entry "+std::to_string(m_entry)+"!");

    cma::DEBUG("EVENT : Set entry for updating ");

    return;
}


void Event::execute(Long64_t entry){
    /* Get the values from the event */
    cma::DEBUG("EVENT : Execute event " );

    // Load data from root tree for this event
    updateEntry(entry);

    // Reset many event-level values
    clear();

    // Get the event weights (for cutflows/histograms)
    initialize_weights();
    cma::DEBUG("EVENT : Setup weights ");

    // Muons
    initialize_muons();
    cma::DEBUG("EVENT : Setup muons ");

    // DNN prediction for each TkMu object
    deepLearningPrediction();

    cma::DEBUG("EVENT : Setup Event ");

    return;
}


void Event::initialize_muons(){
    /* Setup struct of muons and relevant information */
    m_bmtf_muons.clear();
    m_omtf_muons.clear();
    m_emtf_muons.clear();

    // BMTF Muons
    unsigned int nMuons = (*m_mu_bmtf_pt)->size();
    for (unsigned int i=0; i<nMuons; i++){
        Muon mu;
        mu.p4.SetPtEtaPhiM( (*m_mu_bmtf_pt)->at(i),(*m_mu_bmtf_eta)->at(i),(*m_mu_bmtf_phi)->at(i),0);
        mu.charge = (*m_mu_bmtf_charge)->at(i);
        mu.BMTF   = true;
        mu.OMTF   = false;
        mu.EMTF   = false;
        mu.index  = i;
        m_bmtf_muons.push_back(mu);
    }

    // OMTF Muons
    nMuons = (*m_mu_omtf_pt)->size();
    for (unsigned int i=0; i<nMuons; i++){
        Muon mu;
        mu.p4.SetPtEtaPhiM( (*m_mu_omtf_pt)->at(i),(*m_mu_omtf_eta)->at(i),(*m_mu_omtf_phi)->at(i),0);
        mu.charge = (*m_mu_omtf_charge)->at(i);
        mu.BMTF   = false;
        mu.OMTF   = true;
        mu.EMTF   = false;
        mu.index  = i;
        m_omtf_muons.push_back(mu);
    }

    // EMTF Muons
    nMuons = (*m_mu_emtf_pt)->size();
    for (unsigned int i=0; i<nMuons; i++){
        Muon mu;
        mu.p4.SetPtEtaPhiM( (*m_mu_emtf_pt)->at(i),(*m_mu_emtf_eta)->at(i),(*m_mu_emtf_phi)->at(i),0);
        mu.charge = (*m_mu_emtf_charge)->at(i);
        mu.BMTF   = true;
        mu.OMTF   = false;
        mu.EMTF   = false;
        mu.index  = i;
        m_emtf_muons.push_back(mu);
    }

    return;
}


void Event::initialize_tracks(){
    /* Setup vector of tracks */
    m_tracks.clear();

    unsigned int nTracks = (*m_track_pt)->size();
    for (unsigned int i=0; i<nTracks; i++){
        Track tk;

        tk.p4.SetPtEtaPhiM( (*m_track_pt)->at(i), (*m_track_eta)->at(i), (*m_track_phi)->at(i), 0);
        tk.charge = (*m_track_charge)->at(i);
        tk.index  = i;

        tk.d0 = (*m_track_d0)->at(i);
        tk.z0 = (*m_track_z0)->at(i);
        tk.chi2  = (*m_track_chi2)->at(i);
        tk.nstub = (*m_track_nstub)->at(i); 
        tk.fake  = (*m_track_fake)->at(i);
        tk.loose = (*m_track_loose)->at(i);
        tk.genuine = (*m_track_genuine)->at(i);
        tk.unknown = (*m_track_unknown)->at(i);
        tk.combinatoric  = (*m_track_combinatoric)->at(i);
        tk.matchtp_pdgid = (*m_track_matchtp_pdgid)->at(i);
        tk.matchtp_pt  = (*m_track_matchtp_pt)->at(i);
        tk.matchtp_eta = (*m_track_matchtp_eta)->at(i);
        tk.matchtp_phi = (*m_track_matchtp_phi)->at(i);
        tk.matchtp_z0  = (*m_track_matchtp_z0)->at(i);
        tk.matchtp_dxy = (*m_track_matchtp_dxy)->at(i);

        m_tracks.push_back(tk);
    }

    return;
}


void Event::initialize_weights(){
    /* Event weights */
    m_nominal_weight = 1.0;
    return;
}


void Event::deepLearningPrediction(){
    /* Return map of deep learning values */
    cma::DEBUG("EVENT : Calculate DNN features for training ");

    if (!m_DNNtraining && !m_DNNinference){
        cma::ERROR("EVENT : Neither inference nor training -- return");
        return;
    }

    m_tkmuons.clear();

    // loop over tracks
    for (const auto& tk : m_tracks){
        // loop over bmtf muons
        for (const auto& mu : m_bmtf_muons){

            // quality selection
            // -  DeltaR(tk,mu) < 0.8
            // -  pT(track)>2 GeV
            if (tk.p4.DeltaR(mu.p4)>0.8) continue;
            if (tk.p4.Pt()<2.) continue;

            TkMu trkmuon;
            trkmuon.track = tk.index;
            trkmuon.muon  = mu.index;
            if (m_DNNtraining){
                m_deepLearningTool->training(tk,mu);
                trkmuon.features = m_deepLearningTool->features();
            }
            else if (m_DNNinference){
                m_deepLearningTool->inference(tk,mu);
                trkmuon.dnn = m_deepLearningTool->prediction();
            }
            m_tkmuons.push_back( trkmuon );
        } // end loop over muons
    } // end loop over tracks

    return;
}


// -- clean-up
void Event::finalize(){
    /* Delete variables used to access information from TTree */
    cma::DEBUG("EVENT : Finalize -- Clear muons");
    delete m_mu_bmtf_pt;
    delete m_mu_bmtf_eta;
    delete m_mu_bmtf_phi;
    delete m_mu_bmtf_charge;
    delete m_mu_omtf_pt;
    delete m_mu_omtf_eta;
    delete m_mu_omtf_phi;
    delete m_mu_omtf_charge;
    delete m_mu_emtf_pt;
    delete m_mu_emtf_eta;
    delete m_mu_emtf_phi;
    delete m_mu_emtf_charge;

    cma::DEBUG("EVENT : Finalize -- Clear tracks");
    delete m_track_pt;
    delete m_track_eta;
    delete m_track_phi;
    delete m_track_d0;
    delete m_track_z0;
    delete m_track_chi2;
    delete m_track_charge;
    delete m_track_nstub;
    delete m_track_genuine;
    delete m_track_loose;
    delete m_track_unknown;
    delete m_track_fake;
    delete m_track_combinatoric;
    delete m_track_matchtp_pdgid;
    delete m_track_matchtp_pt;
    delete m_track_matchtp_eta;
    delete m_track_matchtp_phi;
    delete m_track_matchtp_z0;
    delete m_track_matchtp_dxy;

    cma::DEBUG("EVENT : Finalize. ");

    return;
}

// THE END
