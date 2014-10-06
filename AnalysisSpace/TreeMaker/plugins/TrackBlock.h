#ifndef __AnalysisSpace_TreeMaker_TrackBlock_h
#define __AnalysisSpace_TreeMaker_TrackBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace vhtm {
  class Track;
}

class TrackBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}

public:
  explicit TrackBlock(const edm::ParameterSet& iConfig);
  virtual ~TrackBlock() {}

  enum {
    kMaxTrack_ = 200
  };

private:
  std::vector<vhtm::Track>* list_;
  int fnTrack_;

  const int verbosity_;
  const edm::InputTag trackTag_;
  const edm::InputTag bsTag_;

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
};
#endif
