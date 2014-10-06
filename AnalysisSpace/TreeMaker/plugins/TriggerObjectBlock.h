#ifndef __AnalysisSpace_TreeMaker_TriggerObjectBlock_h
#define __AnalysisSpace_TreeMaker_TriggerObjectBlock_h

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

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

class TPMERegexp;
namespace vhtm {
  class TriggerObject;
}

class TriggerObjectBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TriggerObjectBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerObjectBlock();

  enum {
    kMaxTriggerObject_ = 100
  };

private:
  int fnTriggerObject_;
  std::vector<vhtm::TriggerObject>* list_;

  const int verbosity_;
  const std::vector<std::string> hltPathsOfInterest_;
  const std::string hltPattern_;
  const double minTrigObjPt_;

  TPMERegexp* re_;
  std::vector<std::string> matchedPathList_;

  const edm::InputTag hltTag_;
  const edm::InputTag objectTag_;

  const edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> objectToken_;
};
#endif
