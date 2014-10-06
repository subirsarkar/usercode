import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
  verbosity = cms.untracked.int32(0),
  patTauSrc = cms.untracked.InputTag('selectedPatTaus'),
  vertexSrc = cms.untracked.InputTag('goodOfflinePrimaryVertices'),
  beamSpotCorr = cms.untracked.bool(True),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot')
)
