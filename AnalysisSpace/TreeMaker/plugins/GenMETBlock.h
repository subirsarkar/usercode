#ifndef __AnalysisSpace_TreeMaker_GenMETBlock_h
#define __AnalysisSpace_TreeMaker_GenMETBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

namespace {
  class GenMET;
}

class GenMETBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit GenMETBlock(const edm::ParameterSet& iConfig);
  virtual ~GenMETBlock() {}

  enum {
    kMaxGenMET_ = 5
  };

private:
  std::vector<vhtm::GenMET>* list_;
  int fnGenMET_;

  const int verbosity_;
  const edm::InputTag genMETTag_;
  const edm::EDGetTokenT<reco::GenMETCollection> genMETToken_;
};
#endif
