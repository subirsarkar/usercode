#ifndef __AnalysisSpace_TreeMaker_VertexBlock_h
#define __AnalysisSpace_TreeMaker_VertexBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace vhtm {
  class Vertex;
}
class VertexBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit VertexBlock(const edm::ParameterSet& iConfig);

  enum {
    kMaxVertex_ = 150
  };

private:
  std::vector<vhtm::Vertex>* list_;
  int fnVertex_;
  int verbosity_;
  edm::InputTag vertexTag_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
};
#endif
