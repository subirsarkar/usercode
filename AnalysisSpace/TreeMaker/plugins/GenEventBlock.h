#ifndef __AnalysisSpace_TreeMaker_GenEventBlock_h
#define __AnalysisSpace_TreeMaker_GenEventBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace {
  class GenEvent;
}
class GenEventBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob() {}

public:
  explicit GenEventBlock(const edm::ParameterSet& iConfig);
  virtual ~GenEventBlock();

private:
  std::vector<vhtm::GenEvent>* list_;
  std::vector<double> *pdfWeights_;

  const int verbosity_;
  const edm::InputTag genEventTag_;
  const bool storePDFWeights_;
  const edm::InputTag pdfWeightsTag_;
  const edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;
  const edm::EDGetTokenT< std::vector<double> > pdfWeightsToken_;
};
#endif
