import FWCore.ParameterSet.Config as cms

trackBlock = cms.EDAnalyzer("TrackBlock",
  verbosity = cms.untracked.int32(0),
  trackSrc = cms.untracked.InputTag('generalTracks'),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot')
)
