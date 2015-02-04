#ifndef __AnalysisSpace_TreeMaker_GenJetBlock_h
#define __AnalysisSpace_TreeMaker_GenJetBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace vhtm {
  class GenJet;
}

class GenJetBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit GenJetBlock(const edm::ParameterSet& iConfig);

  enum {
    kMaxGenJet = 100
  };

private:
  std::vector<vhtm::GenJet>* list_;
  int fnGenJet_;

  const int verbosity_;
  const edm::InputTag genJetTag_;
  const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
};
#endif