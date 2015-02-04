#ifndef __AnalysisSpace_TreeMaker_ElectronBlock_h
#define __AnalysisSpace_TreeMaker_ElectronBlock_h

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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

namespace vhtm {
  class Electron;
}
class ElectronBlock : public edm::EDAnalyzer
{
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

public:
  explicit ElectronBlock(const edm::ParameterSet& iConfig);

  enum {
    kMaxElectron_ = 100
  };

private:
  std::vector<vhtm::Electron>* list_;
  int fnElectron_;

  int verbosity_;
  bool bsCorr_;

  const edm::InputTag bsTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag electronTag_;

  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
};
#endif
