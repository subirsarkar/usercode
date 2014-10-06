#include <memory>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#if 1
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
#endif

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "AnalysisSpace/PattupleMaker/plugins/ElectronUserEmbedded.h"

ElectronsUserEmbedder::ElectronsUserEmbedder(const edm::ParameterSet& iConfig) :
  verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
  electronTag_(iConfig.getParameter<edm::InputTag>("electronTag")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("bsTag", edm::InputTag("offlineBeamSpot"))),
  conversionTag_(iConfig.getUntrackedParameter<edm::InputTag>("conversionTag", edm::InputTag("allConversions"))),
  dcsTag_(iConfig.getUntrackedParameter<edm::InputTag>("dcsTag", edm::InputTag("scalersRawToDigi"))),
  trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag", edm::InputTag("generalTracks"))),
  gsfTrackTag_(iConfig.getUntrackedParameter<edm::InputTag>("gsfTrackTag", edm::InputTag("electronGsfTracks"))),
#if 1
  doMVAPOG_(iConfig.getUntrackedParameter<bool>("doMVAPOG", false)),
#endif
  electronToken_(consumes<pat::ElectronCollection>(electronTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  conversionToken_(consumes<reco::ConversionCollection>(conversionTag_)),
  dcsToken_(consumes<DcsStatusCollection>(dcsTag_)),
  trackToken_(consumes<reco::TrackCollection>(trackTag_)),
  gsfTrackToken_(consumes<reco::GsfTrackCollection>(gsfTrackTag_))
{
#if 1
  edm::FileInPath infileTrig0 = iConfig.getParameter<edm::FileInPath>("infileTrig0");
  edm::FileInPath infileTrig1 = iConfig.getParameter<edm::FileInPath>("infileTrig1");
  edm::FileInPath infileTrig2 = iConfig.getParameter<edm::FileInPath>("infileTrig2");
  edm::FileInPath infileTrig3 = iConfig.getParameter<edm::FileInPath>("infileTrig3");
  edm::FileInPath infileTrig4 = iConfig.getParameter<edm::FileInPath>("infileTrig4");
  edm::FileInPath infileTrig5 = iConfig.getParameter<edm::FileInPath>("infileTrig5");

  edm::FileInPath infileNonTrig0 = iConfig.getParameter<edm::FileInPath>("infileNonTrig0");
  edm::FileInPath infileNonTrig1 = iConfig.getParameter<edm::FileInPath>("infileNonTrig1");
  edm::FileInPath infileNonTrig2 = iConfig.getParameter<edm::FileInPath>("infileNonTrig2");
  edm::FileInPath infileNonTrig3 = iConfig.getParameter<edm::FileInPath>("infileNonTrig3");
  edm::FileInPath infileNonTrig4 = iConfig.getParameter<edm::FileInPath>("infileNonTrig4");
  edm::FileInPath infileNonTrig5 = iConfig.getParameter<edm::FileInPath>("infileNonTrig5");

  if (doMVAPOG_) {
    std::vector<string> wtTrig;
    wtTrig.push_back(infileTrig0.fullPath().data());
    wtTrig.push_back(infileTrig1.fullPath().data());
    wtTrig.push_back(infileTrig2.fullPath().data());
    wtTrig.push_back(infileTrig3.fullPath().data());
    wtTrig.push_back(infileTrig4.fullPath().data());
    wtTrig.push_back(infileTrig5.fullPath().data());
    
    mvaTrig_ = new EGammaMvaEleEstimator();
    mvaTrig_->initialize("BDT",
	                 EGammaMvaEleEstimator::kTrig,
			 true,
			 wtTrig);

    std::vector<string> wtNonTrig;
    wtNonTrig.push_back(infileNonTrig0.fullPath().data());
    wtNonTrig.push_back(infileNonTrig1.fullPath().data());
    wtNonTrig.push_back(infileNonTrig2.fullPath().data());
    wtNonTrig.push_back(infileNonTrig3.fullPath().data());
    wtNonTrig.push_back(infileNonTrig4.fullPath().data());
    wtNonTrig.push_back(infileNonTrig5.fullPath().data());
    
    mvaNonTrig_ = new EGammaMvaEleEstimator();
    mvaNonTrig_->initialize("BDT",
			    EGammaMvaEleEstimator::kNonTrig,
			    true,
			    wtNonTrig);
    
  }
#endif
  produces<pat::ElectronCollection>("");
}
ElectronsUserEmbedder::~ElectronsUserEmbedder()
{
#if 1
  if (doMVAPOG_) {
    delete mvaTrig_;
    delete mvaNonTrig_;
  }
#endif
}
void ElectronsUserEmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::ElectronCollection> electronsHandle;
  iEvent.getByToken(electronToken_, electronsHandle);
  if (!electronsHandle.isValid()) {
    edm::LogError("ElectronUserEmbedded") << "pat::ElectronCollection for " << electronTag_ << " not available!!";
    return;
  }
  const pat::ElectronCollection* electrons = electronsHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);
  if (!vertexHandle.isValid()) {
    edm::LogError("ElectronUserEmbedded") << "reco::VertexCollection for " << vertexTag_ << " not available!!";
    return;
  }
  const reco::VertexCollection* vertexes = vertexHandle.product();

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);
  if (!bsHandle.isValid()) {
    edm::LogError("ElectronUserEmbedded") << "reco::BeamSpot for " << bsTag_ << " not available!!";
    return;
  }
  const reco::BeamSpot &thebs = *bsHandle.product();

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByToken(conversionToken_, hConversions);

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByToken(dcsToken_, dcsHandle);

  float evt_bField = 3.8;
  if (iEvent.isRealData()) {
    if (dcsHandle.isValid()) {
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two-week period
      double currentToBFieldScaleFactor = 2.09237036221512717e-04;
      double current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    }
    else {
      edm::LogError("ElectronUserEmbedded") << "Error >> Failed to get DcsStatusCollection for label: "
                                            << dcsTag_;
    }
  }
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    if (magneticField.isValid()) {
      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    }
    else {
      edm::LogError("ElectronUserEmbedded") << "Error >> Failed to get IdealMagneticFieldRecord";
    }
  }

  // Get the CTF tracks
  edm::Handle<reco::TrackCollection> tracks_h;
  iEvent.getByToken(trackToken_, tracks_h);
  
  // get GSF Tracks
  edm::Handle<reco::GsfTrackCollection> gsftracks_h;
  iEvent.getByToken(gsfTrackToken_, gsftracks_h);

  // Create a new pat collection
  std::auto_ptr<pat::ElectronCollection> electronsUserEmbeddedColl(new pat::ElectronCollection());

  for (auto it = electrons->begin(); it != electrons->end(); ++it) {
    pat::Electron pElectron((*it));
    const reco::GsfElectron* aGsf = static_cast<const reco::GsfElectron*>(&pElectron); 
    if (!aGsf) continue;

    int pfId = 1;

    float nHits = -99;
    if (it->gsfTrack().isNonnull()) {
      const reco::GsfTrackRef tk = it->gsfTrack();
      const reco::HitPattern& p_inner = tk->trackerExpectedHitsInner();
      nHits = p_inner.numberOfHits();
    }
    ConversionFinder convFinder;
    ConversionInfo convInfo = convFinder.getConversionInfo(*aGsf, tracks_h, gsftracks_h, evt_bField);
    double conv_dist = convInfo.dist();
    double conv_dcot = convInfo.dcot();
    int passconv =
      int(!ConversionTools::hasMatchedConversion(*aGsf, hConversions, thebs.position(), true, 2.0, 1e-06, 0));

    float dPhi = pElectron.deltaPhiSuperClusterTrackAtVtx();
    float dEta = pElectron.deltaEtaSuperClusterTrackAtVtx();
    float sihih = pElectron.sigmaIetaIeta();
    float HoE = pElectron.hadronicOverEm();
    if (verbose_) std::cout << "dEta " << dEta 
                            << " dPhi " << dPhi 
                            << " -- dcot " << conv_dcot 
                            << " -- nHits " << nHits 
                            << std::endl;

    pElectron.addUserFloat("nHits", nHits);
    pElectron.addUserFloat("dist", std::abs(conv_dist));
    pElectron.addUserFloat("dcot", std::abs(conv_dcot));
    pElectron.addUserFloat("dPhi", std::abs(dPhi));
    pElectron.addUserFloat("dEta", std::abs(dEta));
    pElectron.addUserFloat("sihih", sihih);
    pElectron.addUserFloat("HoE", HoE);
    pElectron.addUserInt("antiConv", passconv);

    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    if (vertexes->size() > 0) {
      if (pElectron.gsfTrack().isNonnull()) {
        reco::GsfTrackRef tk = it->gsfTrack();
        dxyWrtPV = tk->dxy((*vertexes)[0].position());
        dzWrtPV =  tk->dz((*vertexes)[0].position());
      }
      else if (pElectron.track().isNonnull()) {
        reco::TrackRef tk = it->track();
        dxyWrtPV = tk->dxy((*vertexes)[0].position());
        dzWrtPV =  tk->dz((*vertexes)[0].position());
      }
    }
    pElectron.addUserFloat("dxyWrtPV", dxyWrtPV);
    pElectron.addUserFloat("dzWrtPV", dzWrtPV);

    int myTrigPresel = 0;
    double gsfPt = aGsf->pt();
    if (std::abs(aGsf->superCluster()->eta()) < 1.485) {
      if (aGsf->sigmaIetaIeta() < 0.014 &&
	  aGsf->hadronicOverEm() < 0.15 &&
	  aGsf->dr03TkSumPt()/gsfPt < 0.2 &&
	  aGsf->dr03EcalRecHitSumEt()/gsfPt < 0.2 &&
	  aGsf->dr03HcalTowerSumEt()/gsfPt < 0.2 &&
	  aGsf->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0)
	 myTrigPresel = 1;
    }
    else {
      if (aGsf->sigmaIetaIeta() < 0.035 &&
	  aGsf->hadronicOverEm() < 0.10 &&
	  aGsf->dr03TkSumPt()/gsfPt < 0.2 &&
	  aGsf->dr03EcalRecHitSumEt()/gsfPt < 0.2 &&
	  aGsf->dr03HcalTowerSumEt()/gsfPt < 0.2 &&
	  aGsf->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() == 0)
	myTrigPresel = 1;
    }

    int flagEB = (pElectron.isEB() ? 1 : 0);
    int flagEE = (pElectron.isEE() ? 1 : 0); 
    double pt = pElectron.pt();
  
    int mvaPreselection = passconv && nHits <= 0 && dxyWrtPV < 0.02 && dzWrtPV < 0.1 &&
      ((flagEB 
          && sihih < 0.01
	  && std::abs(dEta) < 0.007
	  && std::abs(dPhi) < 0.15
          && HoE < 0.12
	  && pElectron.dr03TkSumPt()/pt < 0.20
	  && (std::max(pElectron.dr03EcalRecHitSumEt() - 1.0, 0.0))/pt < 0.20
	  && pElectron.dr03HcalTowerSumEt()/pt < 0.20
       ) ||
       (flagEE 
          && sihih < 0.03
	  && std::abs(dEta) < 0.009
	  && std::abs(dPhi) < 0.10
          && HoE < 0.10
	  && pElectron.dr03TkSumPt()/pt < 0.20
	  && (std::max(pElectron.dr03EcalRecHitSumEt() - 1.0, 0.0))/pt < 0.20
	  && pElectron.dr03HcalTowerSumEt()/pt < 0.20
       ));
    pElectron.addUserInt("mvaPreselection", mvaPreselection);
    pElectron.addUserInt("isTriggerElectron", myTrigPresel);

#if 1
    float mva2 = -99;
    float mva3 = -99;
    if (doMVAPOG_) {
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit> > tokenEB 
        = consumes<edm::SortedCollection<EcalRecHit> >(edm::InputTag("reducedEcalRecHitsEB"));
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit> > tokenEE 
        = consumes<edm::SortedCollection<EcalRecHit> >(edm::InputTag("reducedEcalRecHitsEE"));
      EcalClusterLazyTools lazyTools(iEvent, iSetup, tokenEB, tokenEE);

      //EcalClusterLazyTools lazyTools(iEvent, iSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

      edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", hTransientTrackBuilder);
      const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();
  
      mva2 = mvaTrig_->mvaValue(*aGsf, (*vertexes)[0], *transientTrackBuilder, lazyTools, false);
      mva3 = mvaNonTrig_->mvaValue(*aGsf, (*vertexes)[0], *transientTrackBuilder, lazyTools, false);
    }
    pElectron.addUserFloat("mvaPOGTrig", mva2);
    pElectron.addUserFloat("mvaPOGNonTrig", mva3);
#endif
#if 0
    double eta = pElectron.eta();
    double phi = pElectron.phi();

    // iso deposits
    reco::isodeposit::AbsVetos vetosEBPFIdCharged;
    reco::isodeposit::AbsVetos vetosEBPFIdNeutral;
    reco::isodeposit::AbsVetos vetosEBPFIdPhotons;

    reco::isodeposit::AbsVetos vetosEEPFIdCharged;
    reco::isodeposit::AbsVetos vetosEEPFIdNeutral;
    reco::isodeposit::AbsVetos vetosEEPFIdPhotons;

    reco::isodeposit::AbsVetos vetosEBNoPFIdCharged;
    reco::isodeposit::AbsVetos vetosEBNoPFIdNeutral;
    reco::isodeposit::AbsVetos vetosEBNoPFIdPhotons;

    reco::isodeposit::AbsVetos vetosEENoPFIdCharged;
    reco::isodeposit::AbsVetos vetosEENoPFIdNeutral;
    reco::isodeposit::AbsVetos vetosEENoPFIdPhotons;

    // safe recommended: PFId
    vetosEBPFIdCharged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.010));
    vetosEBPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.08));
    vetosEEPFIdCharged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.015));
    vetosEEPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.08));

    // POG recommended: NoPFId
    vetosEENoPFIdCharged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.015));
    vetosEENoPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.08));

    // dr03 
    float chIso03EBPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosEBPFIdCharged).first;
    float nhIso03EBPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosEBPFIdNeutral).first;
    float phIso03EBPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosEBPFIdPhotons).first;
    float nhIsoPU03EBPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosEBPFIdNeutral).first;

    float chIso03EEPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosEEPFIdCharged).first;
    float nhIso03EEPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosEEPFIdNeutral).first;
    float phIso03EEPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosEEPFIdPhotons).first;
    float nhIsoPU03EEPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosEEPFIdNeutral).first;

    float chIso03EBNoPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosEBNoPFIdCharged).first;
    float nhIso03EBNoPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosEBNoPFIdNeutral).first;
    float phIso03EBNoPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosEBNoPFIdPhotons).first;
    float nhIsoPU03EBNoPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosEBNoPFIdNeutral).first;

    float chIso03EENoPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosEENoPFIdCharged).first;
    float nhIso03EENoPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosEENoPFIdNeutral).first;
    float phIso03EENoPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosEENoPFIdPhotons).first;
    float nhIsoPU03EENoPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosEENoPFIdNeutral).first;
    
    // dr04 
    float chIso04EBPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosEBPFIdCharged).first;
    float nhIso04EBPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosEBPFIdNeutral).first;
    float phIso04EBPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosEBPFIdPhotons).first;
    float nhIsoPU04EBPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosEBPFIdNeutral).first;

    float chIso04EEPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosEEPFIdCharged).first;
    float nhIso04EEPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosEEPFIdNeutral).first;
    float phIso04EEPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosEEPFIdPhotons).first;
    float nhIsoPU04EEPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosEEPFIdNeutral).first;

    float chIso04EBNoPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosEBNoPFIdCharged).first;
    float nhIso04EBNoPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosEBNoPFIdNeutral).first;
    float phIso04EBNoPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosEBNoPFIdPhotons).first;
    float nhIsoPU04EBNoPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosEBNoPFIdNeutral).first;

    float chIso04EENoPFId =
      pElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosEENoPFIdCharged).first;
    float nhIso04EENoPFId =
      pElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosEENoPFIdNeutral).first;
    float phIso04EENoPFId =
      pElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosEENoPFIdPhotons).first;
    float nhIsoPU04EENoPFId =
      pElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosEENoPFIdNeutral).first;

    float chIso03 = pfId ? flagEB * chIso03EBPFId   + flagEE * chIso03EEPFId 
                         : flagEB * chIso03EBNoPFId + flagEE * chIso03EENoPFId;
    float nhIso03 = pfId ? flagEB * nhIso03EBPFId   + flagEE * nhIso03EEPFId 
                         : flagEB * nhIso03EBNoPFId + flagEE * nhIso03EENoPFId;
    float phIso03 = pfId ? flagEB * phIso03EBPFId   + flagEE * phIso03EEPFId 
                         : flagEB * phIso03EBNoPFId + flagEE * phIso03EENoPFId;
    float nhIsoPU03 = pfId ? flagEB * nhIsoPU03EBPFId   + flagEE * nhIsoPU03EEPFId 
                           : flagEB * nhIsoPU03EBNoPFId + flagEE * nhIsoPU03EENoPFId;

    float chIso04 = pfId ? flagEB * chIso04EBPFId   + flagEE * chIso04EEPFId 
                         : flagEB * chIso04EBNoPFId + flagEE * chIso04EENoPFId;
    float nhIso04 = pfId ? flagEB * nhIso04EBPFId   + flagEE * nhIso04EEPFId 
                         : flagEB * nhIso04EBNoPFId + flagEE * nhIso04EENoPFId;
    float phIso04 = pfId ? flagEB * phIso04EBPFId   + flagEE * phIso04EEPFId 
                         : flagEB * phIso04EBNoPFId + flagEE * phIso04EENoPFId;
    float nhIsoPU04 = pfId ? flagEB * nhIsoPU04EBPFId   + flagEE * nhIsoPU04EEPFId 
                           : flagEB * nhIsoPU04EBNoPFId + flagEE * nhIsoPU04EENoPFId;

    pElectron.addUserFloat("PFRelIso03",   (chIso03 + nhIso03 + phIso03)/pt);
    pElectron.addUserFloat("PFRelIsoDB03", (chIso03 + std::max(nhIso03 + phIso03 - 0.5*nhIsoPU03, 0.0))/pt);
    pElectron.addUserFloat("PFRelIso04",   (chIso04 + nhIso04 + phIso04)/pt);
    pElectron.addUserFloat("PFRelIsoDB04", (chIso04 + std::max(nhIso04 + phIso04 - 0.5*nhIsoPU04, 0.0))/pt);
   
    std::vector<reco::isodeposit::AbsVetos> list;
    list.push_back(vetosEBPFIdCharged);
    list.push_back(vetosEBPFIdNeutral);
    list.push_back(vetosEBPFIdPhotons);
    list.push_back(vetosEEPFIdCharged);
    list.push_back(vetosEEPFIdNeutral);
    list.push_back(vetosEEPFIdPhotons);
    list.push_back(vetosEBNoPFIdCharged);
    list.push_back(vetosEBNoPFIdNeutral);
    list.push_back(vetosEBNoPFIdPhotons);
    list.push_back(vetosEENoPFIdCharged);
    list.push_back(vetosEENoPFIdNeutral);
    list.push_back(vetosEENoPFIdPhotons);
    for (auto it = list.begin(); it != list.end(); ++it) {
      reco::isodeposit::AbsVetos v = (*it);
      for (auto jt = v.begin(); jt != v.end(); ++jt) {
        delete (*jt);
      }
    }
#endif
    pElectron.addUserInt("pfId", pfId);
    pElectron.addUserFloat("isInRun", iEvent.run());

    electronsUserEmbeddedColl->push_back(pElectron);
  }
  iEvent.put(electronsUserEmbeddedColl);
}
void ElectronsUserEmbedder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
// define this as a plug-in
DEFINE_FWK_MODULE(ElectronsUserEmbedder);
