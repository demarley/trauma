#ifndef EVENTSAVERFLATNTUPLE_H
#define EVENTSAVERFLATNTUPLE_H

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"

#include <memory>
#include <string>
#include <iostream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>

using namespace std;
using namespace edm;


// struct to hold muon data
struct MuonData {
    float pt;
    float eta;
    float phi;
    int q;
};


class EventSaverFlatNtuple : public edm::EDAnalyzer {
public:

  explicit EventSaverFlatNtuple(const edm::ParameterSet& iConfig);
  virtual ~EventSaverFlatNtuple();

  virtual void beginJob() override;
  virtual void endJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  MuonData getMuonData(const l1t::RegionalMuonCand& l1mu) const;

private:
  
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
  int m_process;        // 11/13/211 for single electrons/muons/pions, 6/15 for pions from ttbar/taus, 1 for inclusive
  bool m_debug;         // lots of debug printout statements
  double TP_minPt;      // save TPs with pt > minPt 
  double TP_maxEta;     // save TPs with |eta| < maxEta 
  double TP_maxZ0;      // save TPs with |z0| < maxZ0 
  int m_L1Tk_nparams;   // use 4 or 5 parameter track fit? 
  int m_L1Tk_minNStub;  // require L1 tracks to have >= minNStub (this is mostly for tracklet purposes)
  float m_phi_conv;     // muon local->global phi conversion factor
  double m_mu_maxeta;    // max eta for muons 
  bool m_isMC;

  edm::InputTag L1TrackInputTag;        // L1 track collection
  edm::InputTag MCTruthTrackInputTag;   // MC truth collection
  edm::InputTag MCTruthClusterInputTag;
  edm::InputTag L1StubInputTag;
  edm::InputTag MCTruthStubInputTag;
  edm::InputTag TrackingParticleInputTag;
  edm::InputTag TrackingVertexInputTag;


  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > ttClusterToken_;
  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > ttClusterMCTruthToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;

  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken_;

  edm::InputTag bmtfInputTag;
  edm::InputTag omtfInputTag;
  edm::InputTag emtfInputTag;

  edm::EDGetToken m_bmtfToken;
  edm::EDGetToken m_omtfToken;
  edm::EDGetToken m_emtfToken;

  edm::EDGetTokenT<std::vector<reco::GenParticle>> t_genEvtInfoProd;

  std::vector<unsigned int> m_goodIDs = {6,15,23,24,25};

  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple
  TTree* eventTree;

  std::vector<float> m_trk_pt;
  std::vector<float> m_trk_eta;
  std::vector<float> m_trk_phi;
  std::vector<float> m_trk_d0;           // (filled if L1Tk_nPar==5, else -999)
  std::vector<float> m_trk_z0;
  std::vector<float> m_trk_chi2; 
  std::vector<float> m_trk_stubPtCons;
  std::vector<int>   m_trk_q;
  std::vector<int>   m_trk_nstub;
  std::vector<int>   m_trk_seed;
  std::vector<int>   m_trk_genuine;
  std::vector<int>   m_trk_loose;
  std::vector<int>   m_trk_unknown;
  std::vector<int>   m_trk_combinatoric;
  std::vector<int>   m_trk_fake;         //0 fake; 1 track from PI; 2 secondary track
  std::vector<int>   m_trk_matchtp_pdgid;
  std::vector<float> m_trk_matchtp_pt;
  std::vector<float> m_trk_matchtp_eta;
  std::vector<float> m_trk_matchtp_phi;
  std::vector<float> m_trk_matchtp_z0;
  std::vector<float> m_trk_matchtp_dxy;

  std::vector<float> m_mu_bmtf_pt;
  std::vector<float> m_mu_bmtf_eta;
  std::vector<float> m_mu_bmtf_phi;
  std::vector<int> m_mu_bmtf_q;
  std::vector<float> m_mu_omtf_pt;
  std::vector<float> m_mu_omtf_eta;
  std::vector<float> m_mu_omtf_phi;
  std::vector<int> m_mu_omtf_q;
  std::vector<float> m_mu_emtf_pt;
  std::vector<float> m_mu_emtf_eta;
  std::vector<float> m_mu_emtf_phi;
  std::vector<int> m_mu_emtf_q;

  std::vector<float> m_mc_pt;
  std::vector<float> m_mc_eta;
  std::vector<float> m_mc_phi;
  std::vector<float> m_mc_e;
  std::vector<int> m_mc_pdgId;
  std::vector<int> m_mc_status;
};

#endif
