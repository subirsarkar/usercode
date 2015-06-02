#include <iostream>

#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/MuonBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc", edm::InputTag("selectedPatMuons"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  pfcandTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfCands", edm::InputTag("packedPFCandidates"))),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  muonID_(iConfig.getUntrackedParameter<std::string>("muonID", "GlobalMuonPromptTight")),
  muonToken_(consumes<pat::MuonCollection>(muonTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_)),
  defaultBestMuon_(!iConfig.existsAs<std::string>("customArbitration")),
  bestMuonSelector_(defaultBestMuon_ ? std::string("") : iConfig.getParameter<std::string>("customArbitration"))
{
}
MuonBlock::~MuonBlock() {
  delete list_;
}
void MuonBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Muon>();
  tree->Branch("Muon", "std::vector<vhtm::Muon>", &list_, 32000, 2);
  tree->Branch("nMuon", &fnMuon_, "fnMuon_/I");
}
void MuonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnMuon_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  edm::Handle<pat::MuonCollection> muons;
  bool found = iEvent.getByToken(muonToken_, muons);

  if (found && muons.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(bsToken_, beamSpot);

    unsigned int nMu = muons->size();
    std::vector<bool> gcList(nMu, true);
    for (unsigned int i = 0; i < nMu; ++i) {
      const pat::Muon& mu1 = muons->at(i);
      if (!mu1.track().isNonnull()) {
	gcList[i] = false; 
	continue;
      }
      int nSegments1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
      for (unsigned int j = i+1; j < nMu; ++j) {
        const pat::Muon& mu2 = muons->at(j);
        if (isSameMuon(mu1, mu2)) continue;
        if (!gcList[j] || !mu2.track().isNonnull()) continue;
        int nSegments2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);
        if (nSegments2 == 0 || nSegments1 == 0) continue;
        double sf = muon::sharedSegments(mu1, mu2)/std::min<double>(nSegments1, nSegments2);
	if (sf > 0.499) {
	  if (isBetterMuon(mu1, mu2))
	    gcList[j] = false;
	  else
	    gcList[i] = false;
	}
      }
    }
    
    edm::LogInfo("MuonBlock") << "Total # of Muons: " << muons->size();
    for (unsigned int i = 0; i < nMu; ++i) {
      if (list_->size() == kMaxMuon_) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << list_->size();
	break;
      }
      const pat::Muon& v = muons->at(i);

      vhtm::Muon muon;
      muon.isGlobalMuon  = v.isGlobalMuon() ? true : false;
      muon.isTrackerMuon = v.isTrackerMuon() ? true : false;
      muon.isPFMuon      = v.isPFMuon();
      bool ghostCleaned = gcList[i] || (v.isGlobalMuon() && v.numberOfMatches() >= 2);
      muon.isGhostCleaned = ghostCleaned;

      muon.eta     = v.eta();
      muon.phi     = v.phi();
      muon.pt      = v.pt();
      muon.p       = v.p();
      muon.energy  = v.energy();
      muon.charge  = v.charge();
      muon.passID   = v.muonID(muonID_) ? true : false;
      muon.muonBestTrackType = v.muonBestTrackType();

      double trkd0 = 99;
      double trkdz = 99;
      float normChi2 = 99;

      double dxyWrtPV = 99.;
      double dzWrtPV = 99.;

      double minVtxDist3D = 99.;
      int indexVtx = -1;
      double vertexDistZ = 99.;

      int pixHits = -1;
      int trkHits = -1;
      int muoHits = -1;
      int matches = -1;
      int trackerLayersWithMeasurement = -1;

      reco::TrackRef tk = v.muonBestTrack();
      if (tk.isNonnull()) {
	trkd0 = tk->d0();
	trkdz = tk->dz();
	if (bsCorr_) {
	  if (beamSpot.isValid()) {
	    trkd0 = -(tk->dxy(beamSpot->position()));
	    trkdz = tk->dz(beamSpot->position());
	  }
	  else {
	    edm::LogError("MuonsBlock") << "Error >> Failed to get reco::BeamSpot for label: "
					<< bsTag_;
	  }
	}
	normChi2 = tk->normalizedChi2();
	
	if (primaryVertices.isValid()) {
	  edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();
	  
	  const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
	  dxyWrtPV = tk->dxy(vit.position());
	  dzWrtPV  = tk->dz(vit.position());
	  
	  // Vertex association
	  for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
	    double dxy = tk->dxy(vit->position());
	    double dz = tk->dz(vit->position());
	    double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
	    if (dist3D < minVtxDist3D) {
	      minVtxDist3D = dist3D;
	      indexVtx = static_cast<int>(std::distance(primaryVertices->begin(), vit));
	      vertexDistZ = dz;
	    }
	  }
	}
	else {
	  edm::LogError("MuonBlock") << "Error >> Failed to get reco::VertexCollection for label: "
				     << vertexTag_;
	}
	// Hit pattern
	const reco::HitPattern& hitp = tk->hitPattern();
	pixHits = hitp.numberOfValidPixelHits();
	trkHits = hitp.numberOfValidTrackerHits();
	muoHits = hitp.numberOfValidMuonHits();
	matches = v.numberOfMatches();
	trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();
      }	
      muon.trkD0 = trkd0;
      muon.trkDz = trkdz;
      muon.normChi2 = normChi2;
      muon.dxyPV = dxyWrtPV;
      muon.dzPV  = dzWrtPV;
      muon.vtxDist3D = minVtxDist3D;
      muon.vtxIndex = indexVtx;
      muon.vtxDistZ = vertexDistZ;
      
      muon.pixHits = pixHits;
      muon.trkHits = trkHits;
      muon.muoHits = muoHits;
      muon.matches = matches;
      muon.trackerLayersWithMeasurement = trackerLayersWithMeasurement;
      int numMuonStations = 0;
      unsigned int stationMask = static_cast<unsigned int>(v.stationMask(reco::Muon::SegmentAndTrackArbitration));
      for (int i = 0; i < 8; ++i)  // eight stations, eight bits
	if (stationMask & (1 << i)) ++numMuonStations;
      muon.numMuonStations = numMuonStations;      

      // Isolation
      muon.trkIso   = v.trackIso();
      muon.ecalIso  = v.ecalIso();
      muon.hcalIso  = v.hcalIso();
      muon.hoIso    = v.isolationR03().hoEt;
      
      // PF Isolation
      const reco::MuonPFIsolation& pfIso03 = v.pfIsolationR03();
      muon.sumChargedParticlePtR03 = pfIso03.sumChargedParticlePt;
      muon.sumChargedHadronPtR03   = pfIso03.sumChargedHadronPt;
      muon.sumNeutralHadronEtR03   = pfIso03.sumNeutralHadronEt;
      muon.sumPhotonEtR03 = pfIso03.sumPhotonEt;
      muon.sumPUPtR03 = pfIso03.sumPUPt;
      double absiso = pfIso03.sumChargedHadronPt + std::max(0.0, pfIso03.sumNeutralHadronEt + pfIso03.sumPhotonEt - 0.5 * pfIso03.sumPUPt);
      double iso = absiso/v.p4().pt();
      muon.pfRelIso03 = iso;
      
      const reco::MuonPFIsolation& pfIso04 = v.pfIsolationR04();
      muon.sumChargedParticlePt = pfIso04.sumChargedParticlePt;
      muon.sumChargedHadronPt = pfIso04.sumChargedHadronPt;
      muon.sumNeutralHadronEt = pfIso04.sumNeutralHadronEt;
      muon.sumPhotonEt =  pfIso04.sumPhotonEt;
      muon.sumPUPt = pfIso04.sumPUPt;
      absiso = pfIso04.sumChargedHadronPt + std::max(0.0, pfIso04.sumNeutralHadronEt + pfIso04.sumPhotonEt - 0.5 * pfIso04.sumPUPt);
      iso = absiso/(v.p4().pt());
      muon.pfRelIso04 = iso;
      
      // IP information
      muon.dB = v.dB(pat::Muon::PV2D);
      muon.edB = v.edB(pat::Muon::PV2D);
      
      muon.dB3D = v.dB(pat::Muon::PV3D);
      muon.edB3D = v.edB(pat::Muon::PV3D);
      
      // UW recommendation
      muon.isGlobalMuonPromptTight = muon::isGoodMuon(v, muon::GlobalMuonPromptTight);
      muon.isAllArbitrated         = muon::isGoodMuon(v, muon::AllArbitrated);
      muon.nChambers               = v.numberOfChambers();
      muon.nMatches                = v.numberOfMatches();
      muon.nMatchedStations        = v.numberOfMatchedStations();
      muon.stationMask             = v.stationMask();
      muon.stationGapMaskDistance  = v.stationGapMaskDistance();
      muon.stationGapMaskPull      = v.stationGapMaskPull();
      
      double ptError = tk->ptError()/tk->pt();
      bool muonID = v.isGlobalMuon() && 
	v.isTrackerMuon() && 
	muon.isGlobalMuonPromptTight && 
	muon.isAllArbitrated && 
	std::fabs(dxyWrtPV) < 0.02 && 
	std::fabs(dzWrtPV) < 0.2 && 
	normChi2 < 10 && 
	ptError < 0.1 && 
	muon.trkHits >= 10 && 
	muon.pixHits >= 1 && 
	numMuonStations >= 2 && 
	muon.nMatches >= 1;
      muon.muonID = muonID;
      
      // Vertex information
      const reco::Candidate::Point& vertex = v.vertex();
      muon.vx = vertex.x();
      muon.vy = vertex.y();
      muon.vz = vertex.z();
      
      // Isolation from packed PF candidates 
      std::vector<double> isotemp;
      calcIsoFromPF(v, pfs, 0.15, isotemp);
      muon.isolationMap["c15"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.20, isotemp);
      muon.isolationMap["c20"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.25, isotemp);
      muon.isolationMap["c25"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.30, isotemp);
      muon.isolationMap["c30"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.35, isotemp);
      muon.isolationMap["c35"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.40, isotemp);
      muon.isolationMap["c40"] = isotemp;
      
      isotemp.clear();
      calcIsoFromPF(v, pfs, 0.45, isotemp);
      muon.isolationMap["c45"] = isotemp;
      
      muon.nSegments = v.numberOfMatches(reco::Muon::SegmentArbitration);
      
      list_->push_back(muon);
    }
    fnMuon_ = list_->size();
  }
  else {
    edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon collection for label: "
                               << muonTag_;
  }
}
void MuonBlock::calcIsoFromPF(const pat::Muon& v, 
				const edm::Handle<pat::PackedCandidateCollection>& pfs, 
				double cone, std::vector<double>& iso)
{
  // initialize sums
  double chargedHadSum = 0., 
    chargedSum = 0., 
    neutralSum = 0., 
    photonSum = 0., 
    pileupSum  = 0;

  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0; i < v.numberOfSourceCandidatePtrs(); ++i) 
    footprint.push_back(v.sourceCandidatePtr(i));
  
  // now loop on pf candidates
  for (unsigned int i = 0; i < pfs->size(); ++i) {
    const pat::PackedCandidate& pf = (*pfs)[i];
    int pdgid = std::abs(pf.pdgId());
    double pt = pf.pt();
    if (deltaR(pf, v) < cone) {

      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs, i)) != footprint.end()) continue;
      
      if (pf.charge() == 0) {
        if (pt > 0.5) {
          if (pdgid == 22)
	    photonSum += pt;
          else 
            neutralSum += pt;
        }
      } 
      else if (pf.fromPV() >= 2) {
	chargedSum += pt;
        if (pdgid != 13 && pdgid != 11) chargedHadSum += pt;
      } 
      else
        if (pt > 0.5) pileupSum += pt;
    }
  }
  iso.push_back(chargedHadSum);
  iso.push_back(chargedSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
}
bool MuonBlock::isSameMuon(const pat::Muon& mu1, const pat::Muon& mu2) const {
  return ((&mu1 == &mu2) ||
    (mu1.originalObjectRef() == mu2.originalObjectRef()) ||
    (mu1.reco::Muon::innerTrack().isNonnull() 
     ? mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() 
     : mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack()));
}
bool MuonBlock::isBetterMuon(const pat::Muon &mu1, const pat::Muon &mu2) const {
  if (!defaultBestMuon_) {
    MuonPointerPair pair = {&mu1, &mu2};
    return bestMuonSelector_(pair);
  }
  if (mu2.track().isNull()) return true;
  if (mu1.track().isNull()) return false;
  if (mu1.isPFMuon() != mu2.isPFMuon()) return mu1.isPFMuon();
  if (mu1.isGlobalMuon() != mu2.isGlobalMuon()) return mu1.isGlobalMuon();
  if (mu1.charge() == mu2.charge() && deltaR2(mu1, mu2) < 0.0009) {
    return mu1.track()->ptError()/mu1.track()->pt() < mu2.track()->ptError()/mu2.track()->pt();
  } else {
    int nm1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
    int nm2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);     
    return ((nm1 != nm2) ? nm1 > nm2 : mu1.pt() > mu2.pt());
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
