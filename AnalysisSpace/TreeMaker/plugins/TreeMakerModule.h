#ifndef __AnalysisSpace_TreeMaker_TreeMakerModule_h
#define __AnalysisSpace_TreeMaker_TreeMakerModule_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class TreeMakerModule : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob();

public:
  explicit TreeMakerModule(const edm::ParameterSet& iConfig);
  virtual ~TreeMakerModule() {}

private:
  const int verbosity_;
  const bool createTree_;
};
#endif
