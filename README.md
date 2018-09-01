# Trauma

**Tra**ck-**μ** **Ma**tching with machine learning.

_Given an L1 TT Tracklet and an L1 regional muon candidate, classify the two 
as 'matched' or 'un-matched'._

## Getting Started

Directory | ABOUT
--------- | -----
production/ | Producing flat ntuples from EDM ntuples
training/   | Producing flat ntuples to be used in training (each 'event' represents a track-muon pair)
inference/  | Inference in software with Lightweight NN
fpga/       | Inference in firmware with hls4ml


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

## Questions or comments

Submit an issue or a PR.