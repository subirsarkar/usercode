#ifndef __AnalysisSpace__PattupleMaker__ElectronUserEmbedded_h
#define __AnalysisSpace__PattupleMaker__ElectronUserEmbedded_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class EGammaMvaEleEstimator;
class ElectronsUserEmbedder : public edm::EDProducer {
  public:
    explicit ElectronsUserEmbedder(const edm::ParameterSet&);
    ~ElectronsUserEmbedder();
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
    edm::InputTag electronTag_;
    edm::InputTag vertexTag_;
    edm::InputTag bsTag_;
    edm::InputTag conversionTag_;
    edm::InputTag dcsTag_;
    edm::InputTag trackTag_;
    edm::InputTag gsfTrackTag_;

#if 1
    bool doMVAPOG_;

    EGammaMvaEleEstimator* mvaTrig_;
    EGammaMvaEleEstimator* mvaNonTrig_;
#endif
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<reco::BeamSpot> bsToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionToken_;
    edm::EDGetTokenT<DcsStatusCollection> dcsToken_;
    edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    edm::EDGetTokenT<reco::GsfTrackCollection> gsfTrackToken_;
};
#endif
