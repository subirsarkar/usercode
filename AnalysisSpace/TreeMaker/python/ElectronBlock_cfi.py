import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity = cms.untracked.int32(0),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode  = cms.untracked.bool(False),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('slimmedElectrons'),
  mvaWeightFiles = cms.vstring(
     "EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EB_5_25ns_BDT.weights.xml",
     "EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EE_5_25ns_BDT.weights.xml",
     "EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EB_10_25ns_BDT.weights.xml",
     "EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EE_10_25ns_BDT.weights.xml",
  )
)
