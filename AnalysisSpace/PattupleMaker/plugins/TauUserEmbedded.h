#ifndef __AnalysisSpace__PattupleMaker__TauUserEmbedded_hh
#define __AnalysisSpace__PattupleMaker__TauUserEmbedded_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

class TausUserEmbedded : public edm::EDProducer {
  public:
    explicit TausUserEmbedded(const edm::ParameterSet&);
    ~TausUserEmbedded();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() {}
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() {}
      
    virtual void beginRun(edm::Run&, edm::EventSetup const&) {}
    virtual void endRun(edm::Run&, edm::EventSetup const&) {}
    virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}
    virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

    // ----------member data ---------------------------
    edm::InputTag tauTag_;
    edm::InputTag vertexTag_;

    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
};
#endif
