import FWCore.ParameterSet.Config as cms

from AnalysisSpace.TreeMaker.EventBlock_cfi import eventBlock
from AnalysisSpace.TreeMaker.VertexBlock_cfi import vertexBlock
from AnalysisSpace.TreeMaker.ElectronBlock_cfi import electronBlock
from AnalysisSpace.TreeMaker.MuonBlock_cfi import muonBlock
from AnalysisSpace.TreeMaker.TauBlock_cfi import tauBlock
from AnalysisSpace.TreeMaker.JetBlock_cfi import jetBlock
from AnalysisSpace.TreeMaker.METBlock_cfi import metBlock
from AnalysisSpace.TreeMaker.TriggerBlock_cfi import triggerBlock
from AnalysisSpace.TreeMaker.TriggerObjectBlock_cfi import triggerObjectBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + vertexBlock
 + electronBlock
 + muonBlock
 + tauBlock
 + jetBlock
 + metBlock
 + triggerBlock
 + triggerObjectBlock
)
