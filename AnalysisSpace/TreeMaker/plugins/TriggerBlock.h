#ifndef __AnalysisSpace_TreeMaker_TriggerBlock_h
#define __AnalysisSpace_TreeMaker_TriggerBlock_h

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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class TriggerBlock: public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit TriggerBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerBlock();

private:

  const int verbosity_;

  const edm::InputTag l1Tag_;
  const edm::InputTag hltTag_;
  const edm::InputTag prescaleTag_;
  const std::vector<std::string> hltPathsOfInterest_;

  std::vector<int>* l1physbits_;
  std::vector<int>* l1techbits_;
  std::vector<std::string>* hltpaths_;
  std::vector<int>* hltresults_;
  std::vector<int>* hltprescales_;

  const edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> l1Token_;
  const edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  const edm::EDGetTokenT<pat::PackedTriggerPrescales> prescaleToken_;
};
#endif