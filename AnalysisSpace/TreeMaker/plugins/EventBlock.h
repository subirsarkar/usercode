#ifndef __AnalysisSpace_TreeMaker_EventBlock_h
#define __AnalysisSpace_TreeMaker_EventBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

namespace vhtm {
  class Event;
}
class EventBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob(){}

public:
  explicit EventBlock(const edm::ParameterSet& iConfig);
  virtual ~EventBlock();

private:
  std::vector<vhtm::Event>* list_;
  std::vector<int>* nPU_;
  std::vector<int>* bunchCrossing_;
  std::vector<int>* trueNInt_;

  const int verbosity_;
  const edm::InputTag l1Tag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag trackTag_;
  const edm::InputTag selectedVertexTag_;
  const edm::InputTag puSummaryTag_;
  const edm::InputTag rhoTag_;
  const edm::InputTag rhoNeutralTag_;

  const unsigned int vtxMinNDOF_;
  const double vtxMaxAbsZ_;
  const double vtxMaxd0_;
  const unsigned int numTracks_;
  const double hpTrackThreshold_;

  const edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> l1Token_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::VertexCollection> selectedVertexToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puSummaryToken_;
  const edm::EDGetTokenT<double> rhoToken_;
  const edm::EDGetTokenT<double> rhoNeutralToken_;
};
#endif
