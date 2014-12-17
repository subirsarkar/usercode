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

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/MuonBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  muonTag_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc", edm::InputTag("selectedPatMuons"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  muonID_(iConfig.getUntrackedParameter<std::string>("muonID", "GlobalMuonPromptTight")),
  muonToken_(consumes<pat::MuonCollection>(muonTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_))
{
}
MuonBlock::~MuonBlock() {
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

  edm::Handle<pat::MuonCollection> muons;
  bool found = iEvent.getByToken(muonToken_, muons);

  if (found && muons.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(bsToken_, beamSpot);

    edm::LogInfo("MuonBlock") << "Total # of Muons: " << muons->size();
    for (const pat::Muon& v: *muons) {
      if (list_->size() == kMaxMuon_) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << list_->size();
	break;
      }
      // consider only global muons
      if (!v.isGlobalMuon()) continue;
      reco::TrackRef tk  = v.innerTrack(); // tracker segment only
      reco::TrackRef gtk = v.globalTrack();

      vhtm::Muon muon;
      muon.isTrackerMuon = v.isTrackerMuon() ? true : false;
      muon.isPFMuon      = v.isPFMuon();

      muon.eta     = v.eta();
      muon.phi     = v.phi();
      muon.pt      = v.pt();
      muon.p       = v.p();
      muon.energy  = v.energy();
      muon.charge  = v.charge();

      double trkd0 = tk->d0();
      double trkdz = tk->dz();
      if (bsCorr_) {
        if (beamSpot.isValid()) {
          trkd0 = -(tk->dxy(beamSpot->position()));
          trkdz = tk->dz(beamSpot->position());
        }
        else
          edm::LogError("MuonsBlock") << "Error >> Failed to get reco::BeamSpot for label: "
                                      << bsTag_;
      }
      muon.trkD0      = trkd0;
      muon.trkDz      = trkdz;
      muon.globalChi2 = v.normChi2();
      muon.passID     = v.muonID(muonID_) ? true : false;

      double dxyWrtPV = -99.;
      double dzWrtPV = -99.;
      if (primaryVertices.isValid()) {
        edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
        dxyWrtPV = tk->dxy(vit.position());
        dzWrtPV  = tk->dz(vit.position());
        muon.dxyPV = dxyWrtPV;
        muon.dzPV  = dzWrtPV;

        // Vertex association
        double minVtxDist3D = 9999.;
           int indexVtx = -1;
        double vertexDistZ = 9999.;
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
        muon.vtxDist3D = minVtxDist3D;
        muon.vtxIndex = indexVtx;
        muon.vtxDistZ = vertexDistZ;
      }
      else {
        edm::LogError("MuonBlock") << "Error >> Failed to get reco::VertexCollection for label: "
                                   << vertexTag_;
      }
      // Hit pattern
      const reco::HitPattern& hitp = gtk->hitPattern(); // innerTrack will not provide Muon Hits
      muon.pixHits = hitp.numberOfValidPixelHits();
      muon.trkHits = hitp.numberOfValidTrackerHits();
      muon.muoHits = hitp.numberOfValidMuonHits();
      muon.matches = v.numberOfMatches();
      muon.trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();

      int numMuonStations = 0;
      unsigned int stationMask = static_cast<unsigned int>(v.stationMask(reco::Muon::SegmentAndTrackArbitration));
      for (int i = 0; i < 8; ++i)  // eight stations, eight bits
        if (stationMask & (1<<i)) ++numMuonStations;

      // Isolation
      muon.trkIso   = v.trackIso();
      muon.ecalIso  = v.ecalIso();
      muon.hcalIso  = v.hcalIso();
      muon.hoIso    = v.isolationR03().hoEt;

      // PF Isolation
      const reco::MuonPFIsolation& pfIso03 = v.pfIsolationR03();
      muon.pfChargedIsoR03 = pfIso03.sumChargedParticlePt;
      muon.sumPUPt03 = pfIso03.sumPUPt;
      float absiso = pfIso03.sumChargedParticlePt + std::max(0.0, pfIso03.sumNeutralHadronEt + pfIso03.sumPhotonEt - 0.5 * pfIso03.sumPUPt);
      float iso = absiso/(v.p4().pt());
      muon.pfRelIso03 = iso;

      const reco::MuonPFIsolation& pfIso04 = v.pfIsolationR04();
      muon.pfChargedIsoR04 = pfIso04.sumChargedParticlePt;
      muon.sumPUPt04 = pfIso04.sumPUPt;
      absiso = pfIso04.sumChargedParticlePt + std::max(0.0, pfIso04.sumNeutralHadronEt + pfIso04.sumPhotonEt - 0.5 * pfIso04.sumPUPt);
      iso = absiso/(v.p4().pt());
      muon.pfRelIso04 = iso;

      // IP information
      muon.dB = v.dB(pat::Muon::PV2D);

      // UW recommendation
      muon.isGlobalMuonPromptTight = muon::isGoodMuon(v, muon::GlobalMuonPromptTight);
      muon.isAllArbitrated         = muon::isGoodMuon(v, muon::AllArbitrated);
      muon.nChambers               = v.numberOfChambers();
      muon.nMatches                = v.numberOfMatches();
      muon.nMatchedStations        = v.numberOfMatchedStations();
      muon.stationMask             = v.stationMask();
      muon.stationGapMaskDistance  = v.stationGapMaskDistance();
      muon.stationGapMaskPull      = v.stationGapMaskPull();

      double normalizeChi2 = v.globalTrack()->normalizedChi2();
      double ptError = v.innerTrack()->ptError()/v.innerTrack()->pt();

      bool muonID = v.isGlobalMuon() && 
	v.isTrackerMuon() && 
	muon.isGlobalMuonPromptTight && 
	muon.isAllArbitrated && 
	std::fabs(dxyWrtPV) < 0.02 && 
	std::fabs(dzWrtPV) < 0.2 && 
	normalizeChi2 < 10 && 
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

      list_->push_back(muon);
    }
    fnMuon_ = list_->size();
  }
  else {
    edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon collection for label: "
                               << muonTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
