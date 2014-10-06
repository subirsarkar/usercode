#ifndef __AnalysisSpace_TreeMaker_TauBlock_h
#define __AnalysisSpace_TreeMaker_TauBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

namespace vhtm {
  class Tau;
}
class TauBlock : public edm::EDAnalyzer
{
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

 public:
  explicit TauBlock(const edm::ParameterSet& iConfig);
  virtual ~TauBlock();

  static const reco::PFJetRef& getJetRef(const reco::PFTau& tau);

  enum {
    kMaxTau_ = 100
  };

 private:
  std::vector<vhtm::Tau>* list_;
  int fnTau_;

  const int verbosity_;
  const edm::InputTag tauTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag bsTag_;
  const bool bsCorr_;

  const edm::EDGetTokenT<pat::TauCollection> tauToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
};
#endif
