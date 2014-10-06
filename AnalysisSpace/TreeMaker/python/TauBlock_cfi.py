import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
  verbosity = cms.untracked.int32(0),
  patTauSrc = cms.untracked.InputTag('slimmedTaus'),
  vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  beamSpotCorr = cms.untracked.bool(True),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot')
)
