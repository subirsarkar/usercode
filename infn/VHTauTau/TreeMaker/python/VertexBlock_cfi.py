import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.int32(0),
  vertexSrc = cms.InputTag('goodOfflinePrimaryVertices')
)
