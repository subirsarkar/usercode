import FWCore.ParameterSet.Config as cms

triggerObjectBlock = cms.EDAnalyzer("TriggerObjectBlock",
  verbosity = cms.untracked.int32(0),
  minTrigObjPt = cms.untracked.double(8.0),
  hltInputTag = cms.untracked.InputTag('TriggerResults','','HLT'),
  triggerObjectTag = cms.untracked.InputTag('selectedPatTrigger'),
)
