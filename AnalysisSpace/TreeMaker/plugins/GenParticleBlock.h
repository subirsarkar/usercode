#ifndef __AnalysisSpace_TreeMaker_GenParticleBlock_h
#define __AnalysisSpace_TreeMaker_GenParticleBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

namespace vhtm {
  class GenParticle;
}
class GenParticleBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit GenParticleBlock(const edm::ParameterSet& iConfig);
  virtual ~GenParticleBlock();

  enum {
    kMaxGenParticle_ = 2500
  };

private:
  std::vector<vhtm::GenParticle>* list_;
  int fnGenParticle_;

  const int verbosity_;
  const edm::InputTag genParticleTag_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
};
#endif
