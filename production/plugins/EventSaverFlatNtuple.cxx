/*
Created:        31 August 2018
Last Updated:   31 August 2018

Dan Marley
Texas A&M University
-----

Analyzer to save data to convert EDM to flat ntuple for TRAUMA.
Adapted from
  https://github.com/skinnari/cmssw/blob/TrackletEmulation_937/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc
  - see also: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1TrackletBasedTracking
  - see also: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1Tracklet90X
*/
#include "trauma/production/interface/EventSaverFlatNtuple.h"

using namespace std;
using namespace edm;

EventSaverFlatNtuple::EventSaverFlatNtuple(edm::ParameterSet const& iConfig) :
  config(iConfig){
    m_process        = iConfig.getParameter< int >("process");
    m_debug          = iConfig.getParameter< bool >("debug");
    TP_minPt         = iConfig.getParameter< double >("TP_minPt");
    TP_maxEta        = iConfig.getParameter< double >("TP_maxEta");
    TP_maxZ0         = iConfig.getParameter< double >("TP_maxZ0");
    m_L1Tk_nparams   = iConfig.getParameter< int >("L1Tk_nparams");
    m_L1Tk_minNStub  = iConfig.getParameter< int >("L1Tk_minNStub");
    m_mu_maxeta = iConfig.getParameter< double >("muon_maxeta");
    m_isMC      = iConfig.getParameter<bool>("isMC");           // filling truth branches
    m_phi_conv = 2*M_PI/576.;

    // Muons
    bmtfInputTag = iConfig.getParameter<edm::InputTag>("L1BMTFInputTag");
    omtfInputTag = iConfig.getParameter<edm::InputTag>("L1OMTFInputTag");
    emtfInputTag = iConfig.getParameter<edm::InputTag>("L1EMTFInputTag");

    m_bmtfToken  = consumes< l1t::RegionalMuonCandBxCollection >( bmtfInputTag );
    m_omtfToken  = consumes< l1t::RegionalMuonCandBxCollection >( omtfInputTag );
    m_emtfToken  = consumes< l1t::RegionalMuonCandBxCollection >( emtfInputTag );

    // Tracks
    L1TrackInputTag        = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
    MCTruthTrackInputTag   = iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag");
    MCTruthClusterInputTag = iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag");
    MCTruthStubInputTag    = iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag");

    ttTrackToken_          = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ >>>(L1TrackInputTag);
    ttTrackMCTruthToken_   = consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ >>(MCTruthTrackInputTag);
    ttClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ >>(MCTruthClusterInputTag);
    ttStubMCTruthToken_    = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ >>(MCTruthStubInputTag);

    if (m_isMC)
        t_genEvtInfoProd = consumes<std::vector<reco::GenParticle>>(edm::InputTag("selectedGenParticles", "", "trauma"));
}

EventSaverFlatNtuple::~EventSaverFlatNtuple(){}

void EventSaverFlatNtuple::endJob(){
    /* things to be done at the exit of the event Loop */
    std::cerr << "EventSaverFlatNtuple::endJob" << std::endl;
}

void EventSaverFlatNtuple::beginJob(){
    /* things to be done before entering the event Loop */
    std::cerr << "EventSaverFlatNtuple::beginJob" << std::endl;

    // Create output ntuple
    edm::Service<TFileService> fs;

    eventTree = fs->make<TTree>("training", "training");

    eventTree->Branch("trk_pt",    &m_trk_pt);
    eventTree->Branch("trk_eta",   &m_trk_eta);
    eventTree->Branch("trk_phi",   &m_trk_phi);
    eventTree->Branch("trk_d0",    &m_trk_d0);
    eventTree->Branch("trk_z0",    &m_trk_z0);
    eventTree->Branch("trk_q",     &m_trk_q);
    eventTree->Branch("trk_chi2",  &m_trk_chi2);
    eventTree->Branch("trk_nstub", &m_trk_nstub);
    eventTree->Branch("trk_stubPtCons",   &m_trk_stubPtCons);
    eventTree->Branch("trk_genuine",      &m_trk_genuine);
    eventTree->Branch("trk_loose",        &m_trk_loose);
    eventTree->Branch("trk_unknown",      &m_trk_unknown);
    eventTree->Branch("trk_combinatoric", &m_trk_combinatoric);
    eventTree->Branch("trk_fake",         &m_trk_fake);
    eventTree->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
    eventTree->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
    eventTree->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
    eventTree->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
    eventTree->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
    eventTree->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);

    eventTree->Branch("mu_bmtf_pt",  &m_mu_bmtf_pt );
    eventTree->Branch("mu_bmtf_eta", &m_mu_bmtf_eta);
    eventTree->Branch("mu_bmtf_phi", &m_mu_bmtf_phi);
    eventTree->Branch("mu_bmtf_q",   &m_mu_bmtf_q  );
    eventTree->Branch("mu_omtf_pt",  &m_mu_omtf_pt );
    eventTree->Branch("mu_omtf_eta", &m_mu_omtf_eta);
    eventTree->Branch("mu_omtf_phi", &m_mu_omtf_phi);
    eventTree->Branch("mu_omtf_q",   &m_mu_omtf_q  );
    eventTree->Branch("mu_emtf_pt",  &m_mu_emtf_pt );
    eventTree->Branch("mu_emtf_eta", &m_mu_emtf_eta);
    eventTree->Branch("mu_emtf_phi", &m_mu_emtf_phi);
    eventTree->Branch("mu_emtf_q",   &m_mu_emtf_q  );

    eventTree->Branch("mc_pt",       &m_mc_pt);
    eventTree->Branch("mc_eta",      &m_mc_eta);
    eventTree->Branch("mc_phi",      &m_mc_phi);
    eventTree->Branch("mc_energy",   &m_mc_e);
    eventTree->Branch("mc_id",       &m_mc_pdgId);
    eventTree->Branch("mc_status",   &m_mc_status);

    return;
}


void EventSaverFlatNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    /* Access data from event and store in TTree */
    if (!(m_process==13 || m_process==11 || m_process==211 || m_process==6 || m_process==15 || m_process==1)) {
        std::cout << "The specified m_process is invalid! Exiting..." << std::endl;
        return;
    }

    if ( !(m_L1Tk_nparams==4 || m_L1Tk_nparams==5) ) {
        std::cout << "Invalid number of track parameters, specified m_L1Tk_nparams == " << m_L1Tk_nparams << " but only 4/5 are valid options! Exiting..." << std::endl;
        return;
    }

    // clear variables
    m_trk_pt.clear();
    m_trk_eta.clear();
    m_trk_phi.clear();
    m_trk_d0.clear();
    m_trk_z0.clear();
    m_trk_q.clear();
    m_trk_chi2.clear();
    m_trk_nstub.clear();
    m_trk_stubPtCons.clear();
    m_trk_genuine.clear();
    m_trk_loose.clear();
    m_trk_unknown.clear();
    m_trk_combinatoric.clear();
    m_trk_fake.clear();
    m_trk_matchtp_pdgid.clear();
    m_trk_matchtp_pt.clear();
    m_trk_matchtp_eta.clear();
    m_trk_matchtp_phi.clear();
    m_trk_matchtp_z0.clear();
    m_trk_matchtp_dxy.clear();

    m_mu_bmtf_pt.clear();
    m_mu_bmtf_eta.clear();
    m_mu_bmtf_phi.clear();
    m_mu_bmtf_q.clear();
    m_mu_omtf_pt.clear();
    m_mu_omtf_eta.clear();
    m_mu_omtf_phi.clear();
    m_mu_omtf_q.clear();
    m_mu_emtf_pt.clear();
    m_mu_emtf_eta.clear();
    m_mu_emtf_phi.clear();
    m_mu_emtf_q.clear();

    m_mc_pt.clear();
    m_mc_eta.clear();
    m_mc_phi.clear();
    m_mc_e.clear();
    m_mc_pdgId.clear();
    m_mc_status.clear();


    // Retrieve L1 tracks
    edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
    iEvent.getByToken(ttTrackToken_, TTTrackHandle);

    // Retrieve MC truth - track association maps
    edm::Handle< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ >> MCTruthTTClusterHandle;
    iEvent.getByToken(ttClusterMCTruthToken_, MCTruthTTClusterHandle);
    edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ >> MCTruthTTTrackHandle;
    iEvent.getByToken(ttTrackMCTruthToken_, MCTruthTTTrackHandle);

    // Retrieve Muons
    edm::Handle<BXVector<l1t::RegionalMuonCand>> l1bmtfH;
    edm::Handle<BXVector<l1t::RegionalMuonCand>> l1omtfH;
    edm::Handle<BXVector<l1t::RegionalMuonCand>> l1emtfH;

    iEvent.getByToken(m_bmtfToken, l1bmtfH);
    iEvent.getByToken(m_omtfToken, l1omtfH);
    iEvent.getByToken(m_emtfToken, l1emtfH);

    // Retrieve MC info
    edm::Handle<std::vector<reco::GenParticle>> h_genEvtInfoProd;
    if (m_isMC){
        iEvent.getByToken( t_genEvtInfoProd,h_genEvtInfoProd );
    }

    // ----------------------------------------------------------------------------------------------
    // loop over L1 tracks
    // ----------------------------------------------------------------------------------------------
    if (m_debug) {
        std::cout << std::endl << "Loop over L1 tracks!" << std::endl;
        std::cout << std::endl << "Looking at " << m_L1Tk_nparams << "-parameter tracks!" << std::endl;
    }

    int this_l1track = 0;
    std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
    for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {

        edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);
        this_l1track++;

        float tmp_trk_pt   = iterL1Track->getMomentum(m_L1Tk_nparams).perp();
        float tmp_trk_eta  = iterL1Track->getMomentum(m_L1Tk_nparams).eta();
        float tmp_trk_phi  = iterL1Track->getMomentum(m_L1Tk_nparams).phi();
        float tmp_trk_z0   = iterL1Track->getPOCA(m_L1Tk_nparams).z(); //cm
        float tmp_trk_chi2 = iterL1Track->getChi2(m_L1Tk_nparams);
        int tmp_trk_q      = iterL1Track->getRInv()>0? 1.: -1.;
        float tmp_trk_stubPtCons   = iterL1Track->getStubPtConsistency(m_L1Tk_nparams);
        unsigned int tmp_trk_nstub = iterL1Track->getStubRefs().size();

        float tmp_trk_d0 = -999;
        if (m_L1Tk_nparams == 5) {
            float tmp_trk_x0 = iterL1Track->getPOCA(m_L1Tk_nparams).x();
            float tmp_trk_y0 = iterL1Track->getPOCA(m_L1Tk_nparams).y();
            tmp_trk_d0       = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
        }

        unsigned int tmp_trk_genuine = (MCTruthTTTrackHandle->isGenuine(l1track_ptr));
        unsigned int tmp_trk_loose   = (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr));
        unsigned int tmp_trk_unknown = (MCTruthTTTrackHandle->isUnknown(l1track_ptr));
        unsigned int tmp_trk_combinatoric = (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr));

        if (m_debug) {
            std::cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
                 << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
            if (tmp_trk_genuine) std::cout << " (is genuine)" << std::endl;
            if (tmp_trk_unknown) std::cout << " (is unknown)" << std::endl;
            if (tmp_trk_combinatoric) std::cout << " (is combinatoric)" << std::endl;
        }

        m_trk_q.push_back( tmp_trk_q );
        m_trk_pt .push_back(tmp_trk_pt);
        m_trk_eta.push_back(tmp_trk_eta);
        m_trk_phi.push_back(tmp_trk_phi);
        m_trk_z0 .push_back(tmp_trk_z0);
        m_trk_d0.push_back(tmp_trk_d0);
        m_trk_chi2 .push_back(tmp_trk_chi2);
        m_trk_nstub.push_back(tmp_trk_nstub);
        m_trk_stubPtCons.push_back(tmp_trk_stubPtCons);
        m_trk_genuine.push_back(tmp_trk_genuine);
        m_trk_loose.push_back(tmp_trk_loose);
        m_trk_unknown.push_back(tmp_trk_unknown);
        m_trk_combinatoric.push_back(tmp_trk_combinatoric);


        // ----------------------------------------------------------------------------------------------
        // for studying the fake rate
        // ----------------------------------------------------------------------------------------------
        edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

        int myTP_pdgid = -999;
        float myTP_pt = -999;
        float myTP_eta = -999;
        float myTP_phi = -999;
        float myTP_z0 = -999;
        float myTP_dxy = -999;

        unsigned int myFake = 0;
        if (!my_tp.isNull()) {
            int tmp_eventid = my_tp->eventId().event();

            myFake = (tmp_eventid > 0) ? 2 : 1;

            myTP_pdgid = my_tp->pdgId();
            myTP_pt    = my_tp->p4().pt();
            myTP_eta   = my_tp->p4().eta();
            myTP_phi   = my_tp->p4().phi();
            myTP_z0    = my_tp->vertex().z();

            float myTP_x0 = my_tp->vertex().x();
            float myTP_y0 = my_tp->vertex().y();
            myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

            if (m_debug) {
                  std::cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta()
                       << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
                       << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << std::endl;
            }
        }

        m_trk_fake.push_back(myFake);

        m_trk_matchtp_pdgid.push_back(myTP_pdgid);
        m_trk_matchtp_pt.push_back(myTP_pt);
        m_trk_matchtp_eta.push_back(myTP_eta);
        m_trk_matchtp_phi.push_back(myTP_phi);
        m_trk_matchtp_z0.push_back(myTP_z0);
        m_trk_matchtp_dxy.push_back(myTP_dxy);
    }



    // MUONS
    // Barrel
    for (auto l1mu=l1bmtfH->begin(0); l1mu!=l1bmtfH->end(0);  ++l1mu){ // considering BX = only
        MuonData tmp_data = getMuonData( *l1mu );

        if ( fabs(tmp_data.eta) > m_mu_maxeta) continue;

        m_mu_bmtf_pt.push_back( tmp_data.pt );
        m_mu_bmtf_eta.push_back(tmp_data.eta );
        m_mu_bmtf_phi.push_back(tmp_data.phi );
        m_mu_bmtf_q.push_back(  tmp_data.q );
    }

    // Overlap
    for (auto l1mu=l1omtfH->begin(0); l1mu!=l1omtfH->end(0);  ++l1mu){ // considering BX = only
        MuonData tmp_data = getMuonData( *l1mu );

        if ( fabs(tmp_data.eta) > m_mu_maxeta) continue;

        m_mu_omtf_pt.push_back( tmp_data.pt );
        m_mu_omtf_eta.push_back(tmp_data.eta );
        m_mu_omtf_phi.push_back(tmp_data.phi );
        m_mu_omtf_q.push_back(  tmp_data.q );
    }

    // Endcap
    for (auto l1mu = l1emtfH->begin(0); l1mu != l1emtfH->end(0);  ++l1mu){ // considering BX = only
        MuonData tmp_data = getMuonData( *l1mu );

        if ( fabs(tmp_data.eta) > m_mu_maxeta) continue;

        m_mu_emtf_pt.push_back( tmp_data.pt );
        m_mu_emtf_eta.push_back(tmp_data.eta );
        m_mu_emtf_phi.push_back(tmp_data.phi );
        m_mu_emtf_q.push_back(  tmp_data.q );
    }



    if (m_isMC){
        iEvent.getByToken( t_genEvtInfoProd,h_genEvtInfoProd );

        std::unique_ptr<std::vector<reco::GenParticle> > genCollection( new std::vector<reco::GenParticle> (*h_genEvtInfoProd) );
        for (unsigned int j=0, size=genCollection->size(); j<size; j++){
            reco::GenParticle particle = genCollection->at(j);

            unsigned int absPdgId = std::abs( particle.pdgId() );
            int parent_pdgId(0);
            if (particle.numberOfMothers()>0 && particle.mother(0)!=nullptr)
                parent_pdgId = particle.mother(0)->pdgId();

            // Check that this particle has a PDGID of interest, or that its parent does
            if ( std::find(m_goodIDs.begin(), m_goodIDs.end(), absPdgId) == m_goodIDs.end() &&
                 std::find(m_goodIDs.begin(), m_goodIDs.end(), std::abs(parent_pdgId)) == m_goodIDs.end() )
                continue;


            m_mc_pt.push_back(particle.pt());
            m_mc_eta.push_back(particle.eta());
            m_mc_phi.push_back(particle.phi());
            m_mc_e.push_back(particle.energy());
            m_mc_pdgId.push_back(particle.pdgId());
            m_mc_status.push_back(particle.status());
        } //  end loop over slimmed collection of gen particles
    }


    eventTree->Fill();
} // end of analyze()



MuonData EventSaverFlatNtuple::getMuonData(const l1t::RegionalMuonCand& l1mu) const{
    /* Load the muon data (common to all regions) and return struct */
    MuonData tmp;
    tmp.pt  = l1mu.hwPt() * 0.5;
    tmp.eta = l1mu.hwEta()*0.010875;
    tmp.phi = l1t::MicroGMTConfiguration::calcGlobalPhi( l1mu.hwPhi(),l1mu.trackFinderType(),l1mu.processor() )*m_phi_conv;
    tmp.q   = (l1mu.hwSign()==0) ? 1 : -1;    // rmc.setHwSign(cand->hwSign() == 1 ? 0 : 1 );

    return tmp;
}


DEFINE_FWK_MODULE(EventSaverFlatNtuple);
