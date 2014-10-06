#ifndef __AnalysisSpace_TreeMaker_MuonBlock_h
#define __AnalysisSpace_TreeMaker_MuonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace vhtm {
  class Muon;
}
class MuonBlock : public edm::EDAnalyzer
{
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}

 public:
  explicit MuonBlock(const edm::ParameterSet& iConfig);
  virtual ~MuonBlock();

  enum {
    kMaxMuon_ = 100
  };

 private:
  std::vector<vhtm::Muon>* list_;
  int fnMuon_;

  const int verbosity_;
  const edm::InputTag muonTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag bsTag_;
  const bool bsCorr_;
  const std::string muonID_;

  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
};
#endif
