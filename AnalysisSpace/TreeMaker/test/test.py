import FWCore.ParameterSet.Config as cms
process = cms.Process("TreeMaker")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.infos.threshold = cms.untracked.string("ERROR")
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START71_V1::All'

#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test.root')
)
#--------------------------------------------------
# Analysis Tree Specific
#--------------------------------------------------
process.load("AnalysisSpace.TreeMaker.TreeCreator_cfi")
process.load("AnalysisSpace.TreeMaker.TreeWriter_cfi")
process.load("AnalysisSpace.TreeMaker.TreeContentConfig_cff")

# Primary Vertex Selector
process.selectedPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag('offlineSlimmedPrimaryVertices'),
  cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.0"), 
  filter = cms.bool(True)                                          
)

# Sequence for rhoNeutral
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi import pfAllNeutralHadrons
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import pfAllPhotons
pfNeutralCandPdgIds = []
pfNeutralCandPdgIds.extend(pfAllNeutralHadrons.pdgId.value())
pfNeutralCandPdgIds.extend(pfAllPhotons.pdgId.value())

process.pfNeutralCands = cms.EDFilter("PdgIdPFCandidateSelector",
  src = cms.InputTag('particleFlow'),
  pdgId = cms.vint32(pfNeutralCandPdgIds)
)
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFNeutralJetsForVtxMultReweighting = kt4PFJets.clone(
  src = cms.InputTag('pfNeutralCands'),
  rParam = cms.double(0.6),
  doRhoFastjet = cms.bool(True),
  Rho_EtaMax = cms.double(2.5)
)
process.kt6PFJets = kt4PFJets.clone(
  rParam = cms.double(0.6),
  doRhoFastjet = cms.bool(True),
  Rho_EtaMax = cms.double(2.5)
)

process.p = cms.Path(
  #process.selectedPrimaryVertices *
  #process.kt6PFJets* 
  #process.pfNeutralCands *
  #process.kt6PFNeutralJetsForVtxMultReweighting*
  process.treeCreator*
  process.treeContentSequence*
  process.treeWriter
)

# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
  'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/F6EDDC10-8DFC-E311-BC5D-0025905A60D6.root'
]