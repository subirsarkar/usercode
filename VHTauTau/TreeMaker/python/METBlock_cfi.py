import FWCore.ParameterSet.Config as cms

metBlock = cms.EDAnalyzer("METBlock",
  verbosity = cms.int32(0),
  metSrc    = cms.InputTag("patMETsPF","","PATTuple"),
  corrmetSrc = cms.InputTag("patMETsPFType1Type0","","PATTuple"),
  mvametSrc = cms.InputTag("patMETsPFMVA","","PATTuple")  
)
