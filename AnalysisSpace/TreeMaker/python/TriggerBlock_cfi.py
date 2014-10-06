import FWCore.ParameterSet.Config as cms

triggerBlock = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.untracked.int32(0),
  l1InputTag = cms.untracked.InputTag('gtDigis'),
  hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT'),
  hltPathsOfInterest = cms.vstring ("HLT_DoubleMu",
                                    "HLT_Mu",
                                    "HLT_IsoMu",
                                    "HLT_TripleMu",
                                    "IsoPFTau",
                                    "TrkIsoT",
                                    "HLT_Ele")
)
