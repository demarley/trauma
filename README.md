# Trauma

**Tra**ck-**μ** **M**atching **A**lgorithms

Given an L1 TT Tracklet and an L1 regional muon candidate, classify the two as:
1. Signal (prompt) muon + signal track
2. Signal (prompt) muon + fake track
3. Displaced muon + fake track
4. Fake muon + fake track

## Getting Started

Directory | ABOUT
--------- | -----
production/ | Producing flat ntuples from EDM ntuples
miniAna/    | Producing flat ntuples to be used in training (each 'event' represents a track-muon pair)

_TBA: Instructions for inference in c++ environment (with lwtnn) and hls4ml._

## Resources
This list is not exhaustive, but points to some places I referenced for getting started.

Accessing tracker data:
- CMSSW [L1Tracklets](https://github.com/skinnari/cmssw/blob/TrackletEmulation_937/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker.cc)
- CMSSW Config [L1Tracklets](https://github.com/skinnari/cmssw/blob/TrackletEmulation_937/L1Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker_cfg.py)

Accessing muon data:
- cms-l1t-offline CMSSW [L1TkMuonProducer](https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-integration-CMSSW_10_1_7/L1Trigger/L1TTrackMatch/plugins/L1TkMuonProducer.cc)
- Sven's [MuonConverter](https://gitlab.cern.ch/TrackMuonTriggerCorrelator/Wiki/blob/master/MuonConverter.cc)

Machine Learning:
- [Keras](https://keras.io/)
- [Tensorflow](https://www.tensorflow.org/)
- [hls4ml](https://github.com/hls-fpga-machine-learning/hls4ml)
- [Lightweight NN (lwtnn)](https://github.com/lwtnn/lwtnn)

## Questions or comments

Submit an issue or a PR.