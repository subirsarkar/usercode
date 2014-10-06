#include <memory>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "AnalysisSpace/PattupleMaker/plugins/MuonUserEmbedded.h"

template<typename T>
bool isValidRef(const edm::Ref<T>& ref) {
  return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
}
MuonsUserEmbedded::MuonsUserEmbedded(const edm::ParameterSet& iConfig) :
  verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
  recoMuonTag_(iConfig.getUntrackedParameter<edm::InputTag>("recoMuonTag", edm::InputTag("muons"))),
  pfCandidateTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateTag", edm::InputTag("particleFlow"))),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  muonToken_(consumes<pat::MuonCollection>(muonTag_)),
  recoMuonToken_(consumes<reco::MuonCollection>(recoMuonTag_)),
  pfCandidateToken_(consumes<reco::PFCandidateCollection>(pfCandidateTag_))
{
  produces<pat::MuonCollection>("");
}
MuonsUserEmbedded::~MuonsUserEmbedded() {}
void MuonsUserEmbedded::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Input pat::Muon
  edm::Handle<pat::MuonCollection> muonsHandle;
  iEvent.getByToken(muonToken_, muonsHandle);
  if (!muonsHandle.isValid()) {
    edm::LogError("MuonUserEmbedded") << "pat::MuonCollection for " << muonTag_ << " not available!!";
    return;
  }
  const pat::MuonCollection* muons = muonsHandle.product();

  // reco::Muon
  edm::Handle<reco::MuonCollection> recoMuonsHandle;
  iEvent.getByToken(recoMuonToken_, recoMuonsHandle);
  if (!recoMuonsHandle.isValid()) {
    edm::LogError("MuonUserEmbedded") << "reco::MuonCollection for " << recoMuonTag_ << " not available!!";
    return;
  }
  const reco::MuonCollection* recoMuons = recoMuonsHandle.product();

  // vertices
  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);
  if (!vertexHandle.isValid()) {
    edm::LogError("MuonUserEmbedded") << "reco::VertexCollection for " << vertexTag_ << " not available!!";
    return;
  }
  const reco::VertexCollection* vertexes = vertexHandle.product();

  // ParticleFlow Candidates 
  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByToken(pfCandidateToken_, pfHandle);
  if (!pfHandle.isValid()) {
    edm::LogError("MuonUserEmbedded") << "reco::PFCandidateCollection for " << pfCandidateTag_ << " not available!!";
    return;
  }
  const reco::PFCandidateCollection* pfCandidates = pfHandle.product();

  // create the product object collection
  std::auto_ptr<pat::MuonCollection> muonsUserEmbeddedColl(new pat::MuonCollection());
  for (auto it = muons->begin(); it != muons->end(); ++it) {
    pat::Muon pMuon((*it));

    const reco::Muon* rMuon = 0;
    for (auto jt = recoMuons->begin(); jt != recoMuons->end(); ++jt) {
      if (Geom::deltaR((*jt).p4(), pMuon.p4()) < 1e-03) {
        rMuon = &((*jt));
	edm::LogInfo("MuonUserEmbedded") << "Match to recoMuon" << rMuon->pt();
      }
    }
    if (verbose_) std::cout << "MuonUserEmbedded::recoMuon: " 
			    << rMuon->px() 
                            << ", " << rMuon->py() 
                            << ", " << rMuon->pz() 
                            << std::endl;
    int isPFMuon = 0;
    for (auto jt = pfCandidates->begin(); jt != pfCandidates->end(); ++jt) {
      if ((*jt).particleId() == reco::PFCandidate::mu) {
        reco::MuonRef muonRefToPFMuon = (*jt).muonRef();
        if (muonRefToPFMuon.isNonnull()) {
	  edm::LogInfo("MuonUserEmbedded") << "muonRefToPFMuon: " 
                    << muonRefToPFMuon->px() << ", "
                    << muonRefToPFMuon->py() << ", " 
                    << muonRefToPFMuon->pz(); 
        }
        if (muonRefToPFMuon.isNonnull() &&
            Geom::deltaR(muonRefToPFMuon->p4(), pMuon.p4()) < 1e-04 &&
            (muonRefToPFMuon->isGlobalMuon() || muonRefToPFMuon->isTrackerMuon()))
            isPFMuon = 1;
      }
    }
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    if (vertexes->size() > 0) {
      if (pMuon.isGlobalMuon()) {
        dxyWrtPV = (pMuon.globalTrack())->dxy((*vertexes)[0].position());
        dzWrtPV = (pMuon.globalTrack())->dz((*vertexes)[0].position());
      }
      else if (pMuon.isTrackerMuon()) {
        dxyWrtPV = (pMuon.innerTrack())->dxy((*vertexes)[0].position());
        dzWrtPV = (pMuon.innerTrack())->dz((*vertexes)[0].position());
      }
    }
    pMuon.addUserFloat("dxyWrtPV", dxyWrtPV);
    pMuon.addUserFloat("dzWrtPV", dzWrtPV);
    pMuon.addUserInt("isPFMuon", isPFMuon);

    int globalMuonPromptTight = 0;
    int allArbitrated = 0;
    double normalizeChi2 = 99;
    double ptError = 99;
    int numberOfValidTrackerHits = -99;
    int numberOfValidPixelHits = -99;
    int numMuonStations = -99;
    int numberOfMatches = -99;

    int isGlobalMuon  = pMuon.isGlobalMuon();
    int isTrackerMuon = pMuon.isTrackerMuon();

    bool validRefGlob = isValidRef(pMuon.globalTrack());
    bool validRefInn = isValidRef(pMuon.innerTrack());
    if (validRefGlob && validRefInn) {
      globalMuonPromptTight = muon::isGoodMuon(pMuon, muon::GlobalMuonPromptTight);
      allArbitrated = muon::isGoodMuon(pMuon, muon::AllArbitrated);

      normalizeChi2 = pMuon.globalTrack()->normalizedChi2();
      ptError = (pMuon.innerTrack()->ptError())/(pMuon.innerTrack()->pt());

      const reco::HitPattern& innerTrackHitPattern = pMuon.innerTrack()->hitPattern();
      numberOfValidTrackerHits = innerTrackHitPattern.numberOfValidTrackerHits();
      numberOfValidPixelHits = innerTrackHitPattern.numberOfValidPixelHits();

      int nmst = 0;
      unsigned int stationMask = static_cast<unsigned int>(pMuon.stationMask(reco::Muon::SegmentAndTrackArbitration));
      for (int i = 0; i < 8; ++i)  // eight stations, eight bits
        if (stationMask & (1<<i)) ++nmst;
      numMuonStations = nmst;
      numberOfMatches = pMuon.numberOfMatches();
    }

    int muonID = isGlobalMuon && 
                 isTrackerMuon && 
                 globalMuonPromptTight && 
                 allArbitrated && 
                 (std::abs(dxyWrtPV) < 0.02) && 
                 (std::abs(dzWrtPV) < 0.2) && 
                 (normalizeChi2 < 10) && 
                 (ptError < 0.1) && 
                 (numberOfValidTrackerHits >= 10) && 
                 (numberOfValidPixelHits >= 1) && 
                 (numMuonStations >= 2) && 
                 (numberOfMatches >= 1);

    if (verbose_) std::cout << "isGlobal: " << isGlobalMuon 
                            << " isTracker: " << isTrackerMuon 
                            << " globalMuonPromptTight: " << globalMuonPromptTight
                            << " allArbitrated: " << allArbitrated
                            << " dxyWrtPV: " << dxyWrtPV
                            << " dzWrtPV: " << dzWrtPV 
                            << " normalizeChi2: " << normalizeChi2
                            << " ptError: " << ptError
                            << " numberOfValidTrackerHits: " << numberOfValidTrackerHits
                            << " numberOfValidPixelHits: " << numberOfValidPixelHits
                            << " numMuonStations: " << numMuonStations
                            << " numberOfMatches: " << numberOfMatches
                            << " total ID: " << muonID
                            << std::endl;


    pMuon.addUserInt("isGlobalMuon", isGlobalMuon);
    pMuon.addUserInt("isTrackerMuon", isTrackerMuon);
    pMuon.addUserInt("globalMuonPromptTight", globalMuonPromptTight);
    pMuon.addUserInt("allArbitrated", allArbitrated);
    pMuon.addUserFloat("normalizeChi2", normalizeChi2);
    pMuon.addUserFloat("ptError", ptError);
    pMuon.addUserInt("numberOfValidTrackerHits", numberOfValidTrackerHits);
    pMuon.addUserInt("numberOfValidPixelHits", numberOfValidPixelHits);
    pMuon.addUserInt("numMuonStations", numMuonStations);
    pMuon.addUserInt("numberOfMatches", numberOfMatches);
    pMuon.addUserInt("muonID", muonID);

#if 0    
    double pt = pMuon.pt();
    double eta = pMuon.eta();
    double phi = pMuon.phi();
 
    // iso deposits
    reco::isodeposit::AbsVetos vetosCharged;
    reco::isodeposit::AbsVetos vetosNeutral;
    reco::isodeposit::AbsVetos vetosPhotons;

    vetosCharged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.0001));
    vetosCharged.push_back(new reco::isodeposit::ThresholdVeto(0.0));
    vetosNeutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
    vetosNeutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));
    vetosPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
    vetosPhotons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

    float chIso03 = pMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetosCharged).first;
    float nhIso03 = pMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetosNeutral).first;
    float phIso03 = pMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetosPhotons).first;
    float nhIsoPU03 = pMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosNeutral).first;
    float phIsoPU03 = pMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetosPhotons).first;

    float chIso04 = pMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetosCharged).first;
    float nhIso04 = pMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetosNeutral).first;
    float phIso04 = pMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetosPhotons).first;
    float nhIsoPU04 = pMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosNeutral).first;
    float phIsoPU04 = pMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetosPhotons).first;

    pMuon.addUserFloat("PFRelIso04", (chIso04 + nhIso04 + phIso04)/pt);
    pMuon.addUserFloat("PFRelIso03", (chIso03 + nhIso03 + phIso03)/pt);
    pMuon.addUserFloat("PFRelIsoDB04", 
      (chIso04 + std::max(nhIso04 + phIso04 - 0.5 * 0.5 * (nhIsoPU04 + phIsoPU04), 0.0))/pt);
    pMuon.addUserFloat("PFRelIsoDB03", 
      (chIso03 + std::max(nhIso03 + phIso03 - 0.5 * 0.5 * (nhIsoPU03 + phIsoPU03), 0.0))/pt);

    // cleaning
    std::vector<reco::isodeposit::AbsVetos> list;
    list.push_back(vetosCharged);
    list.push_back(vetosNeutral);
    list.push_back(vetosPhotons);
   
    for (auto it = list.begin(); it != list.end(); ++it) {
      reco::isodeposit::AbsVetos v = (*it);
      for (auto jt = v.begin(); jt != v.end(); ++jt) { 
        delete (*jt);
      }
    }
#endif  
    pMuon.addUserFloat("isInRun", iEvent.run());
    muonsUserEmbeddedColl->push_back(pMuon);
  }
  iEvent.put(muonsUserEmbeddedColl);
}
// ------------ method fills 'descriptions' with the allowed parameters for the module ------------
void MuonsUserEmbedded::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuonsUserEmbedded);
