# define basic process
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TkMu")
 

# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff') ## this needs to match the geometry you are running on
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')     ## this needs to match the geometry you are running on

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))
Source_Files = cms.untracked.vstring(
#  "file:/uscms_data/d3/demarley/correlator/data/singlemuon_noPU_test.root"
   "/store/user/demarley/correlator/singleMuon_PU200.root"
)
process.source = cms.Source("PoolSource", fileNames=Source_Files,
                            inputCommands=cms.untracked.vstring(
                             'keep *',
                             'drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT'
                             )
                            )
process.TFileService = cms.Service("TFileService",fileName = cms.string("output.root"))

# remake stubs 
# ===> IMPORTANT !!! stub window tuning as is by default in CMSSW is incorrect !!! <===
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)


# L1 tracking
process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

## GENERATOR PARTICLES
process.selectedGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
        'drop *',
        'keep status == 3',
        'keep status >= 20 && status <= 100',
        'keep abs(pdgId) == 6 && status >= 20 && status <= 40',
        'keep abs(pdgId) >= 1 && abs(pdgId) <= 5 && status <= 100',
        'keep abs(pdgId) >= 11 && abs(pdgId) <= 18 && status <= 100',
        'keep abs(pdgId) == 23 && status >= 20 && status <= 40',
        'keep abs(pdgId) == 24 && status >= 20 && status <= 100',
        'keep abs(pdgId) == 25 && status >= 20 && status <= 40',
        'keep numberOfMothers() == 1 && abs(mother().pdgId()) == 6 && status >= 20 && status <= 40',
        'keep numberOfMothers() >= 1 && abs(mother().pdgId()) == 24 && status >= 20 && status <= 100',
    )
)

process.tree = cms.EDAnalyzer("EventSaverFlatNtuple",
                                       process = cms.int32(1),
                                       debug   = cms.bool(False),        # printout lots of debug statements
                                       L1Tk_nparams = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.5),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       L1BMTFInputTag  = cms.InputTag("simBmtfDigis","BMTF"),
                                       L1OMTFInputTag  = cms.InputTag("simOmtfDigis","OMTF"),
                                       L1EMTFInputTag  = cms.InputTag("simEmtfDigis","EMTF"),
                                       muon_maxeta = cms.double(5),      # only save TPs with |eta| < X
                                       isMC = cms.bool(False),
               )



process.FEVToutput_step = cms.EndPath(process.tree)

process.schedule = cms.Schedule(process.TTClusterStub,process.TTTracksWithTruth,process.FEVToutput_step)

