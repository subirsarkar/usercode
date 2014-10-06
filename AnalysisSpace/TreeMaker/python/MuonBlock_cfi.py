import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  verbosity = cms.untracked.int32(0),
  muonSrc = cms.untracked.InputTag('selectedPatMuons'),
  vertexSrc = cms.untracked.InputTag('goodOfflinePrimaryVertices'),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  beamSpotCorr = cms.untracked.bool(True),
  muonID = cms.untracked.string('GlobalMuonPromptTight')
)
