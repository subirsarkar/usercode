import FWCore.ParameterSet.Config as cms

process = cms.Process("PATN")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 1

# import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# to run in un-scheduled mode uncomment the following lines
process.options.allowUnscheduled = cms.untracked.bool(True)

# parse arguments
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register ('isMC',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Sample Type: MC or data")

options.register ('channel',
                  'ett',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Desired channel: mtt, ett or none")    

options.register ('includeSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include Sim. Default: False")

options.register ('includePatTrig',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Include PatTrig. Default: False")

import sys
print sys.argv

if len(sys.argv) > 0:
  last = sys.argv.pop()
  sys.argv.extend(last.split(","))
  print sys.argv

if hasattr(sys, "argv") == True:
  options.parseArguments()
isMC = options.isMC
channel = options.channel
includeSim = options.includeSim
includePatTrig = options.includePatTrig
print 'Using channel: %s' % channel

# Customize PAT object creation
postfix = ""
jetAlgo = "AK5"
excludeFromTopProjection=['Tau']
from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,
          runPF2PAT=True, 
          jetAlgo=jetAlgo, 
          runOnMC=isMC, 
          postfix=postfix,
          excludeFromTopProjection=excludeFromTopProjection, 
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=True)

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)
process.muonVariables = cms.EDProducer('MuonsUserEmbedded',
    muonTag = cms.InputTag("selectedPatMuons"),
    vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS")
)
process.electronVariables = cms.EDProducer('ElectronsUserEmbedder',
    electronTag = cms.InputTag("selectedPatElectrons"),
    vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS")
)
process.tauVariables = cms.EDProducer('TausUserEmbedded',
       tauTag = cms.InputTag("selectedPatTaus"),
    vertexTag = cms.InputTag("offlinePrimaryVerticesWithBS")
)
## Add new patME collecions

## Embed tracks
# PAT Muons
process.patMuons.embedTrack = cms.bool(True) # Embed tracks for muon ID cuts to be done offline

# PAT Electrons
process.patElectrons.embedTrack = cms.bool(True)
process.patElectrons.embedGsfTrack = cms.bool(True)

# PAT Taus
process.patTaus.embedLeadPFCand = cms.bool(True)
process.patTaus.embedSignalPFCands = cms.bool(True)
process.patTaus.embedIsolationPFCands = cms.bool(True)
process.patTaus.embedLeadTrack = cms.bool(True)
process.patTaus.embedSignalTracks = cms.bool(True)
process.patTaus.embedIsolationTracks = cms.bool(True)
process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(True)
process.patTaus.embedIsolationPFNeutralHadrCands = cms.bool(True)
process.patTaus.embedIsolationPFGammaCands = cms.bool(True)
process.patTaus.embedSignalPFChargedHadrCands = cms.bool(True)
process.patTaus.embedSignalPFNeutralHadrCands = cms.bool(True)
process.patTaus.embedSignalPFGammaCands = cms.bool(True)
process.patTaus.embedLeadPFChargedHadrCand = cms.bool(True)
process.patTaus.embedLeadPFNeutralCand = cms.bool(True)

# good Vertex collection
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True)
)

# Trigger and Trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process)

# skim
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True)
)

simpleCutsWP95 = "(userFloat('nHits')<=1"+ \
                 " && (" + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.80 && "+ \
                 "          userFloat('dEta') <0.007 && userFloat('HoE') <0.15)"   + \
                 " || "  + \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.70 && "+ \
                 "          userFloat('dEta') <0.010 && userFloat('HoE') <0.07)"   + \
                 "     )"+ \
                 ")"

process.skimmedMuons = cms.EDFilter("PATMuonSelector",
  src = cms.InputTag("muonVariables"),
  cut = cms.string("pt >= 20. && abs(eta) < 2.1 && isGlobalMuon && (trackIso + ecalIso +hcalIso)/pt < 0.3"),
  filter = cms.bool(True)
)

process.skimmedElectrons = cms.EDFilter("PATElectronSelector",
  src = cms.InputTag("electronVariables"),
  cut = cms.string("pt >= 20. && abs(eta) < 2.1 && ( " + simpleCutsWP95 + " ) && (trackIso + ecalIso +hcalIso)/pt < 0.3"),
  filter = cms.bool(True)
)

process.skimmedTaus = cms.EDFilter("PATTauSelector",
  src = cms.InputTag("tauVariables"),
  cut = cms.string("pt >= 20. && abs(eta) < 2.3 && tauID('decayModeFinding') > 0.5"),
  filter = cms.bool(True)
)

process.numTaus = cms.EDFilter("PATCandViewCountFilter",
  src = cms.InputTag("skimmedTaus"),
  maxNumber = cms.uint32(2000),
  minNumber = cms.uint32(2),
  filter = cms.bool(True)
)

if (channel == "mtt"):
  process.theSkim = cms.Sequence(
		process.goodOfflinePrimaryVertices +
    		process.skimmedMuons + 
		process.skimmedTaus +
		process.numTaus)
elif (channel == "ett"):
  process.theSkim = cms.Sequence(
		process.goodOfflinePrimaryVertices +
    		process.skimmedElectrons + 
		process.skimmedTaus +
		process.numTaus)
else:
  process.theSkim = cms.Sequence()

# Event content
#from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent,patTriggerEventContent,patTriggerL1RefsEventContent
#process.out.outputCommands.extend(patExtraAodEventContent)
#process.out.outputCommands.extend(patTriggerEventContent)
#process.out.outputCommands.extend(patTriggerL1RefsEventContent)
process.out.outputCommands = [
    'drop *',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_genParticles_*_*',
    'keep *_TriggerResults_*_HLT',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep *_gtDigis_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_TriggerResults_*_PAT',
    'keep *_patTriggerEvent_*_*',
    'keep patMuons_muonVariables_*_*',
    'keep patElectrons_electronVariables_*_*',
    'keep patTaus_tauVariables_*_*',
    'keep *_selectedPatJets_*_*',
    'keep *_patMETs_*_*',
    'keep *_selectedPatPFParticles_*_*',
#    'keep *_selectedPatPhotons_*_*',
    'keep *_patTrigger_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep PileupSummaryInfos_*_*_*',
    'keep recoVertexs_offlinePrimaryVertices_*_*',
    'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
    'keep recoVertexs_goodOfflinePrimaryVertices_*_*'
]

## let it run
process.p = cms.Path(
    process.scrapingVeto +
    process.goodOfflinePrimaryVertices +
    process.patDefaultSequence+
    process.muonVariables +
    process.electronVariables +
    process.tauVariables
    #+process.theSkim
)

## remove MC matching from the default sequence
if not isMC:
  removeMCMatching(process, ['All'])
  removeMCMatching(process, ['METs'], "TC")
  removeMCMatching(process, ['METs'], "PF")
  process.patDefaultSequence.remove(process.patJetPartonMatch)
  #process.patDefaultSequence.remove(process.patJetPartonMatchAK5PF)
  #process.patDefaultSequence.remove(process.patJetGenJetMatchAK5PF)
  process.patDefaultSequence.remove(process.patJetFlavourId)
  process.patDefaultSequence.remove(process.patJetPartons)
  process.patDefaultSequence.remove(process.patJetPartonAssociation)
  #process.patDefaultSequence.remove(process.patJetPartonAssociationAK5PF)
  process.patDefaultSequence.remove(process.patJetFlavourAssociation)
  #process.patDefaultSequence.remove(process.patJetFlavourAssociationAK5PF)
  runOnData(process)

if not isMC:
  process.patTriggerEvent.processName = '*'
  if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'

from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
process.source.fileNames = filesRelValProdTTbarAODSIM
#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(
#  '/store/mc/Fall13dr/WH_ZH_HToTauTau_M-125_13TeV_pythia6/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/087890E7-6A80-E311-8AA4-001E67398110.root'
#  ),
#                            skipEvents = cms.untracked.uint32(0)
#)
process.out.fileName = cms.untracked.string('patTuple.root')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')
