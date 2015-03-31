import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity = cms.untracked.int32(0),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode  = cms.untracked.bool(False),
  offlineBeamSpot = cms.untracked.InputTag('offlineBeamSpot'),
  vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
  electronSrc = cms.untracked.InputTag('slimmedElectrons'),
  mvaWeightFiles = cms.vstring(
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml",
        "EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml"
  )
)
