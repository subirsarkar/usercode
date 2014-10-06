#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TPRegexp.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "AnalysisSpace/TreeMaker/plugins/TriggerBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

static const unsigned int NmaxL1AlgoBit = 128;
static const unsigned int NmaxL1TechBit = 64;

// Constructor
TriggerBlock::TriggerBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  l1Tag_(iConfig.getUntrackedParameter<edm::InputTag>("l1InputTag", edm::InputTag("gtDigis"))),
  hltTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults","","HLT"))),
  prescaleTag_(iConfig.getUntrackedParameter<edm::InputTag>("prescaleInputTag", edm::InputTag("patTrigger"))),
  hltPathsOfInterest_(iConfig.getParameter<std::vector<std::string> >("hltPathsOfInterest")),
  l1Token_(consumes<L1GlobalTriggerReadoutRecord>(l1Tag_)),
  hltToken_(consumes<edm::TriggerResults>(hltTag_)),
  prescaleToken_(consumes<pat::PackedTriggerPrescales>(prescaleTag_))
{
}
TriggerBlock::~TriggerBlock() {
  delete l1physbits_;
  delete l1techbits_;
  delete hltpaths_;
  delete hltresults_;
  delete hltprescales_;
}
void TriggerBlock::beginJob()
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);

  l1physbits_ = new std::vector<int>();
  tree->Branch("l1physbits", "vector<int>", &l1physbits_);

  l1techbits_ = new std::vector<int>();
  tree->Branch("l1techbits", "vector<int>", &l1techbits_);

  hltpaths_ = new std::vector<std::string>();
  tree->Branch("hltpaths", "vector<string>", &hltpaths_);

  hltresults_ = new std::vector<int>();
  tree->Branch("hltresults", "vector<int>", &hltresults_);

  hltprescales_ = new std::vector<int>();
  tree->Branch("hltprescales", "vector<int>", &hltprescales_);
}
void TriggerBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vectors
  l1physbits_->clear();
  l1techbits_->clear();
  hltpaths_->clear();
  hltresults_->clear();
  hltprescales_->clear();

  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  bool found = iEvent.getByToken(l1Token_, l1GtReadoutRecord);
  if (found && l1GtReadoutRecord.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: "
                                 << l1Tag_;

    for (unsigned int i = 0; i < NmaxL1AlgoBit; ++i) 
      l1physbits_->push_back(l1GtReadoutRecord->decisionWord()[i] ? 1 : 0);

    for (unsigned int i = 0; i < NmaxL1TechBit; ++i) 
      l1techbits_->push_back(l1GtReadoutRecord->technicalTriggerWord()[i] ? 1 : 0 );
  }
  else 
    edm::LogError("TriggerBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label: "
                                  << l1Tag_;

  edm::Handle<edm::TriggerResults> triggerResults;
  found = iEvent.getByToken(hltToken_, triggerResults);
  if (found && triggerResults.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained " << hltTag_;

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(prescaleToken_, triggerPrescales);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
    for (unsigned int i = 0; i < triggerResults->size(); ++i) {
      std::string path = names.triggerName(i);
      if (hltPathsOfInterest_.size()) {
        int nmatch = 0;
        for (auto kt: hltPathsOfInterest_) {
          nmatch += TPRegexp(kt).Match(path);
        }
        if (!nmatch) continue;
      }
      hltpaths_->push_back(path);
      hltprescales_->push_back(triggerPrescales->getPrescaleForIndex(i));
      hltresults_->push_back((triggerResults->accept(i) ? 1 : 0));
    }
  } 
  else {
    edm::LogError("TriggerBlock") << "Failed to get TriggerResults for label: "
                                  << hltTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerBlock);
