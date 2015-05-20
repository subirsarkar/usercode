#ifndef __AnalysisSpace_TreeMaker_PackedCandidateBlock_h
#define __AnalysisSpace_TreeMaker_PackedCandidateBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

namespace vhtm {
  class PackedPFCandidate;
}

class PackedPFCandidateBlock : public edm::EDAnalyzer {
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}
  void calcIsoFromPF(const pat::PackedCandidate& v, edm::Handle<pat::PackedCandidateCollection>& pfs, double cone, std::vector<double>& iso);
 public:
  explicit PackedPFCandidateBlock(const edm::ParameterSet& iConfig);
  virtual ~PackedPFCandidateBlock() {}
  
  enum {
    kMaxPackedPFCandidate = 250
  };
 private:
  std::vector<vhtm::PackedPFCandidate>* list_;
  int fnPackedPFCandidate_;
  
  int verbosity_;
  const edm::InputTag pfcandTag_;
  std::vector<int> pdgTosave_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};
#endif
