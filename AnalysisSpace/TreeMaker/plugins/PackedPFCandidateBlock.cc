#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "AnalysisSpace/TreeMaker/plugins/PackedPFCandidateBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
PackedPFCandidateBlock::PackedPFCandidateBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfcandTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCands", edm::InputTag("packedPFCandidates"))), 
  pdgTosave_(iConfig.getParameter<std::vector<int>>("pdgTosave")),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_))
{}
void PackedPFCandidateBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::PackedPFCandidate>();
  tree->Branch("PackedPFCandidate", "std::vector<vhtm::PackedPFCandidate>", &list_, 32000, 2);
  tree->Branch("nPackedPFCandidate", &fnPackedPFCandidate_, "fnPackedPFCandidate_/I");
}
void PackedPFCandidateBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnPackedPFCandidate_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  
  if (pfs.isValid()) {
    edm::LogInfo("PackedPFCandidateBlock") << "Total # PackedPFCandidate: " << pfs->size();
    for (pat::PackedCandidate const& v: *pfs) {
      if (list_->size() == kMaxPackedPFCandidate) {
      	edm::LogInfo("PackedPFCandidateBlock") << "Too many PackedPFCandidates, fnPackedPFCandidate = " 
					       << fnPackedPFCandidate_; 
      	break;
      }
      
      int pdg = std::abs(v.pdgId());
      if (std::find(pdgTosave_.begin(), pdgTosave_.end(), pdg) == pdgTosave_.end() || v.pt() <= 2.) continue;
      
      vhtm::PackedPFCandidate pfCand;
      
      pfCand.pt = v.pt();
      pfCand.eta = v.eta();
      pfCand.phi = v.phi();
      pfCand.energy = v.energy();
      
      pfCand.pdgId = v.pdgId();
      pfCand.charge = v.charge();
      
      pfCand.vx = v.vx();
      pfCand.vz = v.vz();
      pfCand.vz = v.vz();
      
      pfCand.fromPV = v.fromPV();
      pfCand.dxy = v.dxy();
      pfCand.dz = v.dz();
      pfCand.dxyError = v.dxyError();   
      pfCand.dzError = v.dzError();   
      
      std::vector<double> isotemp;   
      calcIsoFromPF(0.30, pfs, v, isotemp);
      pfCand.isolationMap["c30"] = isotemp;
      
      list_->push_back(pfCand);
    }
    fnPackedPFCandidate_ = list_->size();
  }
  else {
    edm::LogError("PackedPFCandidateBlock") << "Error >> Failed to get pat::PackedPFCandidate for label: " 
					    << pfcandTag_;
  }
}
void PackedPFCandidateBlock::calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, 
					   const pat::PackedCandidate& v, std::vector<double>& iso)
{
  // initialize sums
  double chargedHadSum = 0., 
    chargedParticleSum = 0., 
    neutralSum = 0., 
    photonSum = 0., 
    pileupSum  = 0;

  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0; i < v.numberOfSourceCandidatePtrs(); ++i)
    footprint.push_back(v.sourceCandidatePtr(i));
 
  // now loop on pf candidates
  for (unsigned int i = 0; i < pfs->size(); ++i) {
    const pat::PackedCandidate& pf = (*pfs)[i];
    double dRcone = deltaR(v, pf);
    if (dRcone < cone) {
      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs, i)) != footprint.end()) continue;

      double pt = pf.pt(); 
      int pdg = std::abs(pf.pdgId());
      if (pf.charge() == 0) {
        if (pt > 0.5 && dRcone > 0.01) {
          if (pdg == 22)
            photonSum += pt;
          else 
            neutralSum += pt;
        }
      } 
      else {
        if (pt > 0.2 && dRcone > 0.0001) { 
          if (pf.fromPV() >= 2) {
	    chargedParticleSum += pt;
            if (pdg != 13 && pdg != 11) chargedHadSum += pt;
          } 
	  else 
            pileupSum += pt;
        }
      }
    }
  }
  if (verbosity_) std::cout << "isoValues: (" << chargedHadSum << "," 
			    << neutralSum << "," << photonSum << "," 
			    << pileupSum << ")" 
			    << std::endl;
  iso.push_back(chargedHadSum);
  iso.push_back(chargedParticleSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PackedPFCandidateBlock);
