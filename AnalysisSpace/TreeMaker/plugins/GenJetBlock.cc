#include <iostream>

#include "TTree.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisSpace/TreeMaker/plugins/GenJetBlock.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

GenJetBlock::GenJetBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genJetTag_(iConfig.getUntrackedParameter<edm::InputTag>("genJetSrc", edm::InputTag("ak5GenJets"))),
  genJetToken_(consumes<reco::GenJetCollection>(genJetTag_))   
{
}
void GenJetBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::GenJet>();
  tree->Branch("GenJet", "std::vector<vhtm::GenJet>", &list_, 32000, 2);
  tree->Branch("nGenJet", &fnGenJet_, "fnGenJet_/I");
}
void GenJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  list_->clear();
  fnGenJet_ = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenJetCollection> genJets;
    bool found = iEvent.getByToken(genJetToken_, genJets);

    if (found && genJets.isValid()) {
      edm::LogInfo("GenJetBlock") << "Total # GenJets: " << genJets->size();
      for (auto it = genJets->begin(); it != genJets->end(); ++it) {
        if (fnGenJet_ == kMaxGenJet) {
	  edm::LogInfo("GenJetBlock") << "Too many Gen Jets, fnGenJet = " << fnGenJet_;
	  break;
        }
	vhtm::GenJet jet;

        // fill in all the vectors
        jet.eta = it->eta();
        jet.phi = it->phi();
        jet.p   = it->p();
        jet.pt  = it->pt();
        double energy = it->energy();
        jet.energy = energy;
        jet.emf    = (energy > 0) ? it->emEnergy()/energy : 0;
        jet.hadf   = (energy > 0) ? it->hadEnergy()/energy : 0;
 
        list_->push_back(jet);
      }
      fnGenJet_ = list_->size();
    }
    else {
      edm::LogError("GenJetBlock") << "Error >> Failed to get GenJetCollection for label: "
                                   << genJetTag_;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
