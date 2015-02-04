#ifndef __AnalysisSpace_TreeMaker_METBlock_h
#define __AnalysisSpace_TreeMaker_METBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/MET.h"

namespace vhtm {
  class MET;
}
class METBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit METBlock(const edm::ParameterSet& iConfig);
  virtual ~METBlock();

  enum {
    kMaxMET_ = 5
  };

  void fillMET(const edm::Event& iEvent,
               const edm::EventSetup& iSetup,
               std::vector<vhtm::MET>* list,
               int& fnMET,
               const edm::InputTag& iTag,  
               const edm::EDGetTokenT<pat::METCollection>& token);

private:
  std::vector<vhtm::MET>* pfList_;
  int fnPFMET_;
#if 0
  std::vector<vhtm::MET>* corrList_;
  int fnCorrMET_;

  std::vector<vhtm::MET>* mvaList_;
  int fnMVAMET_;
#endif
  const int verbosity_;
  const edm::InputTag pfMETTag_;
#if 0
  const edm::InputTag corrMETTag_;
  const edm::InputTag mvaMETTag_;
#endif
  const edm::EDGetTokenT<pat::METCollection> pfMETToken_;
#if 0
  const edm::EDGetTokenT<pat::METCollection> corrMETToken_;
  const edm::EDGetTokenT<pat::METCollection> mvaMETToken_;
#endif
};
#endif
