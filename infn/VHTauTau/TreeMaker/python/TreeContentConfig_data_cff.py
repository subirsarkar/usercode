import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.EventBlock_cfi import eventBlock 
from VHTauTau.TreeMaker.VertexBlock_cfi import vertexBlock 
from VHTauTau.TreeMaker.ElectronBlock_cfi import electronBlock 
from VHTauTau.TreeMaker.MuonBlock_cfi import muonBlock 
from VHTauTau.TreeMaker.TauBlock_cfi import tauBlock 
from VHTauTau.TreeMaker.JetBlock_cfi import jetBlock 
from VHTauTau.TreeMaker.METBlock_cfi import metBlock 
from VHTauTau.TreeMaker.TriggerBlock_cfi import triggerBlock 
from VHTauTau.TreeMaker.TriggerObjectBlock_cfi import triggerObjectBlock

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
