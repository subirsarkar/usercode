#ifndef __AnalysisSpace__PattupleMaker__MuonUserEmbedded_h
#define __AnalysisSpace__PattupleMaker__MuonUserEmbedded_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

class MuonsUserEmbedded : public edm::EDProducer {
  public:
    explicit MuonsUserEmbedded(const edm::ParameterSet&);
    ~MuonsUserEmbedded();

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
    bool verbose_;
    edm::InputTag muonTag_;
    edm::InputTag vertexTag_;
    edm::InputTag recoMuonTag_;
    edm::InputTag pfCandidateTag_;

    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<reco::MuonCollection> recoMuonToken_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateToken_;
};
#endif
