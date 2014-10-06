## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# verbose flags for the PF2PAT modules
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service("Tracer")

runOnMC = True

if runOnMC:
    from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
    process.source.fileNames = filesRelValProdTTbarAODSIM
else:
    from PhysicsTools.PatAlgos.patInputFiles_cff import filesSingleMuRECO
    process.source.fileNames = filesSingleMuRECO
    process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )

# load the PAT config (for the PAT only part)
#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
addMETCollection(process, labelName='patMETsTypeIcorrected', metSource='pfType1CorrectedMet')
addMETCollection(process, labelName='patMETTC', metSource='tcMet')

# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
postfix = ""
jetAlgo = "AK5"
#Define Objects to be excluded from Top Projection. Default is Tau, so objects are not cleaned for taus
excludeFromTopProjection=['Tau']

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix,excludeFromTopProjection=excludeFromTopProjection)
# to run second PF2PAT+PAT with different postfix uncomment the following lines
# and add the corresponding sequence to path
#postfix2 = "PFlow2"
#jetAlgo2="AK7"
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo2, runOnMC=True, postfix=postfix2)

# to use tau-cleaned jet collection uncomment the following:
#getattr(process,"pfNoTau"+postfix).enable = True

if not runOnMC:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatchingPF2PAT( process, '' )

# Add PF2PAT output to the created file
from PhysicsTools.PatAlgos.patEventContent_cff import \
patEventContentNoCleaning,patExtraAodEventContent,patTriggerEventContent,patTriggerL1RefsEventContent,patEventContentTriggerMatch
process.out.outputCommands.extend(patExtraAodEventContent)
process.out.outputCommands.extend(patTriggerEventContent)
process.out.outputCommands.extend(patTriggerL1RefsEventContent)
process.out.outputCommands.extend(patEventContentTriggerMatch)
process.out.outputCommands.extend(['drop recoGenJets_*_*_*'])
process.out.outputCommands.extend(['keep recoVertexs_goodOfflinePrimaryVertices_*_*',
                               'keep recoPFCandidates_particleFlow_*_*',
                               'keep PileupSummaryInfos_*_*_*',
                               'keep recoVertexs_goodOfflinePrimaryVertices_*_*',
                               'keep *_offlineBeamSpot_*_*',
                               'keep L1GlobalTriggerReadoutRecord_*_*_*',
                               'keep *_ak5GenJets_*_*',
                               'keep *_kt6GenJets_*_*',
                               'keep *_genMetTrue_*_*',
                               #'keep *_kt6PFJets_*_*',
                               'keep recoPFJets_ak5PFJets_*_*',
                               'keep *_*GsfElectrons*_*_*',
                               'keep *_electronGsfTracks_*_*'])
process.p = cms.Path(
  process.patDefaultSequence
)

## ------------------------------------------------------
# In addition you usually want to change the following
# parameters:
## ------------------------------------------------------
#
# process.GlobalTag.globaltag = ... ## (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
# ##
# process.source.fileNames = ... ## (e.g. 'file:AOD.root')
#
process.maxEvents.input = 10
# ##
# process.out.outputCommands = [ ... ] ## (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
# ##
process.out.fileName = 'patTuple_PATandPF2PAT.root'
# ##
# process.options.wantSummary = False ## (to suppress the long output at the end of the job)

#process.prune(verbose=True)
