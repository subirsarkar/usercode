## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.options.allowUnscheduled = cms.untracked.bool(True)

## to run in un-scheduled mode uncomment the following lines
#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

#from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
#addMETCollection(process, labelName='patMETsTypeIcorrected', metSource='pfType1CorrectedMet')
#addMETCollection(process, labelName='patMETTC', metSource='tcMet')

runOnMC=True
# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
postfix = ""
jetAlgo = "AK5"
#Define Objects to be excluded from Top Projection. Default is Tau, so objects are not cleaned for taus
excludeFromTopProjection=['Tau']

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT(process, runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, excludeFromTopProjection=excludeFromTopProjection)

#from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
#from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

## uncomment the following lines to add ak5PFJetsCHS to your PAT output
#postfixAK5PFCHS = 'Copy'
#addJetCollection(
#   process,
#   postfix = postfixAK5PFCHS,
#   labelName = 'AK5PFCHS',
#   jetSource = cms.InputTag('ak5PFJetsCHS'),
#   jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2')
#)
#process.out.outputCommands.append( 'drop *_selectedPatJetsAK5PFCHS%s_caloTowers_*'%( postfixAK5PFCHS ) )

# uncomment the following lines to add ak5PFJets to your PAT output
#addJetCollection(
#   process,
#   labelName = 'AK5PF',
#   jetSource = cms.InputTag('ak5PFJets'),
#   jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'),
#   btagDiscriminators = [
#       'jetBProbabilityBJetTags'
#     , 'jetProbabilityBJetTags'
#     , 'trackCountingHighPurBJetTags'
#     , 'trackCountingHighEffBJetTags'
#     , 'simpleSecondaryVertexHighEffBJetTags'
#     , 'simpleSecondaryVertexHighPurBJetTags'
#     , 'combinedSecondaryVertexBJetTags'
#   ],
#)
#process.out.outputCommands.append( 'drop *_selectedPatJetsAK5PF_caloTowers_*' )
process.patJets.addJetID=True
process.patJets.jetIDMap="ak5JetID"
#process.patJets.useLegacyJetMCFlavour=True

#from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
#process.source.fileNames = filesRelValProdTTbarAODSIM
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
  '/store/mc/Fall13dr/WH_ZH_HToTauTau_M-125_13TeV_pythia6/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/087890E7-6A80-E311-8AA4-001E67398110.root'
  ),
                            skipEvents = cms.untracked.uint32(0)
)
process.maxEvents.input = 10

# Output
from PhysicsTools.PatAlgos.patEventContent_cff import patExtraAodEventContent
process.out.fileName = cms.untracked.string('patTuple_mc.root') 
#process.out.outputCommands.extend(patExtraAodEventContent)
#process.out.outputCommands += ['drop recoGenJets_*_*_*' ]
#process.out.outputCommands += [
#                              'keep recoVertexs_goodOfflinePrimaryVertices_*_*',
#                               'keep PileupSummaryInfos_*_*_*',
#                               'keep *_ak5GenJets_*_*',
#                               'keep *_kt6GenJets_*_*',
#                               'keep *_genMetTrue_*_*',
#                               #'keep *_kt6PFJets_*_*',
#                               'keep recoPFJets_ak5PFJets_*_*',
#                               'keep *_*GsfElectrons*_*_*',
#                               'keep *_electronGsfTracks_*_*'
#]
print process.out.outputCommands.dumpPython()
process.load('PhysicsTools.PatAlgos.patSequences_cff')
process.p = cms.Path(
#  process.kt6PFJets*
  process.patDefaultSequence
)
process.options.wantSummary = False
