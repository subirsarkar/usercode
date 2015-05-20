#ifndef __AnalysisSpace_TreeMaker_MuonBlock_h
#define __AnalysisSpace_TreeMaker_MuonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace vhtm {
  class Muon;
}
class MuonBlock : public edm::EDAnalyzer
{
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

 public:
  explicit MuonBlock(const edm::ParameterSet& iConfig);
  virtual ~MuonBlock();
  void calcIsoFromPF(const pat::Muon& v, const edm::Handle<pat::PackedCandidateCollection>& pfs, double cone, std::vector<double>& iso);     

  enum {
    kMaxMuon_ = 100
  };

 private:
  std::vector<vhtm::Muon>* list_;
  int fnMuon_;

  const int verbosity_;
  const bool keepOnlyGlobalMuons_;

  const edm::InputTag muonTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag bsTag_;
  const edm::InputTag pfcandTag_;
  const bool bsCorr_;
  const std::string muonID_;

  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};
#endif
