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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

namespace vhtm {
  class TriggerObject;
}

class TriggerObjectBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);

public:
  explicit TriggerObjectBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerObjectBlock();
  static void printObjectInfo(const pat::TriggerObjectStandAlone& obj); 

  enum {
    kMaxTriggerObject_ = 200
  };

private:
  int fnTriggerObject_;
  std::vector<vhtm::TriggerObject>* list_;

  const int verbosity_;
  const double minTrigObjPt_;
  HLTConfigProvider hltConfig_;

  const edm::InputTag hltTag_;
  const edm::InputTag objectTag_;

  const edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> objectToken_;
};
#endif
