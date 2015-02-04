#include <iostream>
#include <algorithm>

#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/JetBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

// Constructor
JetBlock::JetBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  jetTag_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc",  edm::InputTag("selectedPatJets"))),
  jetToken_(consumes<pat::JetCollection>(jetTag_))
{}
void JetBlock::beginJob()
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::Jet>();
  tree->Branch("Jet", "std::vector<vhtm::Jet>", &list_, 32000, 2);
  tree->Branch("nJet", &fnJet_, "fnJet_/I");
}
void JetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnJet_ = 0;

  edm::Handle<pat::JetCollection> jets;
  bool found = iEvent.getByToken(jetToken_, jets);

  if (found && jets.isValid()) {
    unsigned int njets = jets->size();
    edm::LogInfo("JetBlock") << "Total # PAT Jets: " << njets;
    for (const pat::Jet& jet: *jets) {
      if (list_->size() == kMaxJet_) {
        edm::LogInfo("JetBlock") << "Too many PAT Jets, fnJet = " << list_->size();
        break;
      }
      retpf.set(false);
      int passjetLoose = (pfjetIDLoose(jet, retpf)) ? 1 : 0;

      retpf.set(false);
      int passjetTight = (pfjetIDTight(jet, retpf)) ? 1 : 0;

      vhtm::Jet jobj;

      // fill in all the vectors
      jobj.eta    = jet.eta();
      jobj.phi    = jet.phi();
      jobj.pt     = jet.pt();
      jobj.pt_raw = jet.correctedJet("Uncorrected").pt();
      jobj.energy = jet.energy();
      jobj.energy_raw    = jet.correctedJet("Uncorrected").energy();
      jobj.partonFlavour = jet.partonFlavour();

      jobj.chargedEmEnergyFraction     = jet.chargedEmEnergyFraction();
      jobj.chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
      jobj.chargedMuEnergyFraction     = jet.chargedMuEnergyFraction();
      jobj.electronEnergyFraction      = jet.electronEnergy()/jet.energy();
      jobj.muonEnergyFraction          = jet.muonEnergyFraction();
      jobj.neutralEmEnergyFraction     = jet.neutralEmEnergyFraction();
      jobj.neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
      jobj.photonEnergyFraction        = jet.photonEnergyFraction();

      jobj.chargedHadronMultiplicity  = jet.chargedHadronMultiplicity();
      jobj.chargedMultiplicity        = jet.chargedMultiplicity();
      jobj.electronMultiplicity       = jet.electronMultiplicity();
      jobj.muonMultiplicity           = jet.muonMultiplicity();
      jobj.neutralHadronMultiplicity  = jet.neutralHadronMultiplicity();
      jobj.neutralMultiplicity        = jet.neutralMultiplicity();
      jobj.photonMultiplicity         = jet.photonMultiplicity();

      jobj.nConstituents = jet.numberOfDaughters();

      jobj.combinedSecondaryVertexBTag = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
      jobj.combinedInclusiveSecondaryVertexBTag 
                                       = jet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
      jobj.passLooseID = passjetLoose;
      jobj.passTightID = passjetTight;

      for (const std::pair<std::string, float>& pa: jet.getPairDiscri())
	jobj.discrimap[pa.first] = pa.second;

      list_->push_back(jobj);
    }
    fnJet_ = list_->size();
  }
  else {
    edm::LogError("JetBlock") << "Error >> Failed to get pat::Jet collection for label: "
                              << jetTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetBlock);
