import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.untracked.int32(1),
  vertexSrc = cms.untracked.InputTag('offlineSlimmedPrimaryVertices')
)
