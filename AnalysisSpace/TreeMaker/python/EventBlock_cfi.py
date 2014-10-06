import FWCore.ParameterSet.Config as cms

eventBlock = cms.EDAnalyzer("EventBlock",
  verbosity = cms.untracked.int32(0),
  l1InputTag = cms.untracked.InputTag('gtDigis'),
  vertexTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  vertexMinimumNDOF = cms.untracked.uint32(4),
  vertexMaxAbsZ = cms.untracked.double(24.),
  vertexMaxd0 = cms.untracked.double(2.),
  #trkInputTag = cms.untracked.InputTag('generalTracks'),
  #numTracks = cms.untracked.uint32(10),
  hpTrackThreshold = cms.untracked.double(0.25),
  puSummaryInputTag = cms.untracked.InputTag('addPileupInfo'),
  #selectedVtxInputTag = cms.untracked.InputTag('selectedPrimaryVertices'),
  selectedVtxInputTag = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  rhoInputTag = cms.untracked.InputTag('kt6PFJets','rho'),                         
  rhoNeutralInputTag = cms.untracked.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho')                         
)
