import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  verbosity = cms.untracked.int32(0),
  muonSrc = cms.untracked.InputTag('slimmedMuons'),
  vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  beamSpotCorr = cms.untracked.bool(True),
  muonID = cms.untracked.string('GlobalMuonPromptTight')
)
