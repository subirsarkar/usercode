#ifndef __VHTauTau_TreeMaker_JetBlock_h
#define __VHTauTau_TreeMaker_JetBlock_h

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
#include "DataFormats/PatCandidates/interface/Jet.h"

namespace vhtm {
  class Jet;
}
class JetBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit JetBlock(const edm::ParameterSet& iConfig);

  enum {
    kMaxJet_ = 100
  };
private:
  std::vector<vhtm::Jet>* list_;
  int fnJet_;

  const int verbosity_;
  const edm::InputTag jetTag_;
  const edm::EDGetTokenT<pat::JetCollection> jetToken_;
};
#endif
