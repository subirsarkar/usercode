#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/GenParticleBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

using namespace std;

GenParticleBlock::GenParticleBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genParticleTag_(iConfig.getUntrackedParameter<edm::InputTag>("genParticleSrc", edm::InputTag("genParticles"))),
  genParticleToken_(consumes<reco::GenParticleCollection>(genParticleTag_))
{
}
GenParticleBlock::~GenParticleBlock() {
}
void GenParticleBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::GenParticle>();
  tree->Branch("GenParticle", "std::vector<vhtm::GenParticle>", &list_, 32000, 2);
  tree->Branch("nGenParticle", &fnGenParticle_, "fnGenParticle_/I");
}
void GenParticleBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnGenParticle_ = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    bool found = iEvent.getByToken(genParticleToken_, genParticles);
    if (found && genParticles.isValid()) {
      edm::LogInfo("GenParticleBlock") << "Total # GenParticles: " << genParticles->size();
      if (verbosity_ > 0) {
	std::cout << setprecision(2);
        std::cout << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
		  << std::endl;
      }
      for (auto it = genParticles->begin(); it != genParticles->end(); ++it) {
        if (fnGenParticle_ == kMaxGenParticle_) {
           edm::LogInfo("GenParticleBlock") << "Too many GenParticles, fnGenParticle = "
                                            << fnGenParticle_;
           break;
        }
        int pdgid = it->pdgId();
        double pt = it->pt();

        vhtm::GenParticle gp;

        // fill in all the vectors
        gp.eta       = it->eta();
        gp.phi       = it->phi();
        gp.p         = it->p();
        gp.px        = it->px();
        gp.py        = it->py();
        gp.pz        = it->pz();
        gp.pt        = pt;
        gp.energy    = it->energy();
        gp.pdgId     = pdgid;
        gp.vx        = it->vx();
        gp.vy        = it->vy();
        gp.vz        = it->vz();
        gp.status    = it->status();
        gp.charge    = it->charge();
        gp.numDaught = it->numberOfDaughters();
        gp.numMother = it->numberOfMothers();

        // First mother
        int idx = -1;
        for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
          if (it->mother() == &(*mit)) {
            idx = std::distance(genParticles->begin(), mit);
            break;
          }
        }
        gp.motherIndex = idx;

        // Mothers
        ostringstream mID;
        for (size_t j = 0; j < it->numberOfMothers(); ++j) {
          const reco::Candidate* m = it->mother(j);
          for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
            if (m == &(*mit) ) {
              int idx = std::distance(genParticles->begin(), mit);
              gp.motherIndices.push_back(idx);
 	      mID << " " << idx; 
              break;
            }
          }
        }

        // Daughters
        ostringstream dID;
        for (size_t j = 0; j < it->numberOfDaughters(); ++j) {
          const reco::Candidate* d = it->daughter(j);
          for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
            if (d == &(*mit) ) {
              int idx = std::distance(genParticles->begin(), mit);
              gp.daughtIndices.push_back(idx);
	      dID << " " << idx; 
              break;
            }
          }
        }

        // Dump
        if (verbosity_ > 0) {
          string ms = mID.str();
	  if (!ms.length()) ms = " -";

          string ds = dID.str();
	  if (!ds.length()) ds = " -";

	  std::cout << setw(4)  << fnGenParticle_++
		    << setw(8)  << it->status()
		    << setw(10) << it->pdgId()
                    << setw(8)  << it->charge()
		    << setw(10) << it->eta()
		    << setw(9)  << it->phi()
		    << setw(9)  << it->pt()
		    << setw(9)  << it->energy()
		    << setw(16) << ms 
		    << ds
		    << std::endl;
        }
        // add particle to the list
        list_->push_back(gp);
      }
      fnGenParticle_ = list_->size();
    }
    else {
      edm::LogError("GenParticleBlock") << "Error >> Failed to get GenParticleCollection for label: "
                                        << genParticleTag_;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleBlock);
