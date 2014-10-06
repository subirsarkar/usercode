import FWCore.ParameterSet.Config as cms

from AnalysisSpace.TreeMaker.EventBlock_cfi import eventBlock
#from AnalysisSpace.TreeMaker.GenEventBlock_cfi import genEventBlock
from AnalysisSpace.TreeMaker.VertexBlock_cfi import vertexBlock
#from AnalysisSpace.TreeMaker.GenJetBlock_cfi import genJetBlock
from AnalysisSpace.TreeMaker.JetBlock_cfi import jetBlock
from AnalysisSpace.TreeMaker.ElectronBlock_cfi import electronBlock
from AnalysisSpace.TreeMaker.METBlock_cfi import metBlock
#from AnalysisSpace.TreeMaker.GenMETBlock_cfi import genMETBlock
from AnalysisSpace.TreeMaker.MuonBlock_cfi import muonBlock
from AnalysisSpace.TreeMaker.TauBlock_cfi import tauBlock
from AnalysisSpace.TreeMaker.GenParticleBlock_cfi import genParticleBlock
from AnalysisSpace.TreeMaker.TriggerBlock_cfi import triggerBlock
from AnalysisSpace.TreeMaker.TriggerObjectBlock_cfi import triggerObjectBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + triggerBlock
 + triggerObjectBlock
 + vertexBlock
# + genEventBlock
 + genParticleBlock
# + genJetBlock
# + genMETBlock
 + electronBlock
  + muonBlock
 + tauBlock
 + metBlock
 + jetBlock
)
