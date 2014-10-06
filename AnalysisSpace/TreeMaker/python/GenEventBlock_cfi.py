import FWCore.ParameterSet.Config as cms

genEventBlock = cms.EDAnalyzer("GenEventBlock",
  verbosity = cms.untracked.int32(0),
  GenEventInfoInputTag = cms.untracked.InputTag('generator'),
  PDFWeightsInputTag = cms.untracked.InputTag('pdfWeights','cteq66'),
  StorePDFWeights = cms.untracked.bool(False)
)
