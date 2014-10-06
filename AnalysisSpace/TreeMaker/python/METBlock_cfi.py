import FWCore.ParameterSet.Config as cms

metBlock = cms.EDAnalyzer('METBlock',
  verbosity = cms.untracked.int32(0),
  metSrc = cms.untracked.InputTag('slimmedMETs'),
  corrmetSrc = cms.untracked.InputTag('patMETsTypeIcorrected'),
  mvametSrc = cms.untracked.InputTag('patPFMetByMVA')
)
