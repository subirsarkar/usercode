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

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig) :
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
    for (auto it = muons->begin(); it != muons->end(); ++it) {
      if (fnMuon_ == kMaxMuon_) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << fnMuon_;
	break;
      }
      // consider only global muons
      if (!it->isGlobalMuon()) continue;
      reco::TrackRef tk  = it->innerTrack(); // tracker segment only
      reco::TrackRef gtk = it->globalTrack();

      vhtm::Muon muon;
      muon.isTrackerMuon = (it->isTrackerMuon()) ? true : false;
      muon.isPFMuon      = it->userInt("isPFMuon") ? true : false;

      muon.eta     = it->eta();
      muon.phi     = it->phi();
      muon.pt      = it->pt();
      muon.ptError = tk->ptError();
      muon.p       = it->p();
      muon.energy  = it->energy();
      muon.charge  = it->charge();

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
      muon.trkD0Error = tk->d0Error();
      muon.trkDz      = trkdz;
      muon.trkDzError = tk->dzError();
      muon.globalChi2 = it->normChi2();
      muon.passID     = (it->muonID(muonID_)) ? true : false;

      if (primaryVertices.isValid()) {
        edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        auto vit = primaryVertices->begin(); // Highest sumPt vertex
        muon.dxyPV = tk->dxy(vit->position());
        muon.dzPV  = tk->dz(vit->position());

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
            indexVtx = int(std::distance(primaryVertices->begin(), vit));
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
      muon.matches = it->numberOfMatches();
      muon.trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();

      // Isolation
      muon.trkIso   = it->trackIso();
      muon.ecalIso  = it->ecalIso();
      muon.hcalIso  = it->hcalIso();
      muon.hoIso    = it->isolationR03().hoEt;
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();
      muon.relIso = reliso;

      // PF Isolation
      muon.pfChargedIsoR03 = it->pfIsolationR03().sumChargedParticlePt;
      muon.pfChargedIsoR04 = it->pfIsolationR04().sumChargedParticlePt;

#if 0
      muon.pfRelIso = it->userFloat("PFRelIsoDB04v2");
#endif
      // IP information
      muon.dB    = it->dB(pat::Muon::PV2D);
      muon.edB   = it->edB(pat::Muon::PV2D);
      muon.dB3d  = it->dB(pat::Muon::PV3D);
      muon.edB3d = it->edB(pat::Muon::PV3D);

      // UW recommendation
      muon.isGlobalMuonPromptTight = muon::isGoodMuon(*it, muon::GlobalMuonPromptTight);
      muon.isAllArbitrated         = muon::isGoodMuon(*it, muon::AllArbitrated);
      muon.nChambers               = it->numberOfChambers();
      muon.nMatches                = it->numberOfMatches();
      muon.nMatchedStations        = it->numberOfMatchedStations();
      muon.stationMask             = it->stationMask();
      muon.stationGapMaskDistance  = it->stationGapMaskDistance();
      muon.stationGapMaskPull      = it->stationGapMaskPull();
      muon.muonID                  = it->userInt("muonID");

      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      muon.vx = vertex.x();
      muon.vy = vertex.y();
      muon.vz = vertex.z();

#if 0
      muon.idMVA          = it->userFloat("muonIdMVA");
      muon.isoRingsMVA    = it->userFloat("muonIsoRingsMVA");
      muon.isoRingsRadMVA = it->userFloat("muonIsoRingsRadMVA");
      muon.idIsoCombMVA   = it->userFloat("muonIdIsoCombMVA");

      muon.pfRelIso03v1   = it->userFloat("PFRelIso03v1");
      muon.pfRelIso03v2   = it->userFloat("PFRelIso03v2");
      muon.pfRelIsoDB03v1 = it->userFloat("PFRelIsoDB03v1");
      muon.pfRelIsoDB03v2 = it->userFloat("PFRelIsoDB03v2");

      muon.pfRelIso04v1   = it->userFloat("PFRelIso04v1");
      muon.pfRelIso04v2   = it->userFloat("PFRelIso04v2");
      muon.pfRelIsoDB04v1 = it->userFloat("PFRelIsoDB04v1");
      muon.pfRelIsoDB04v2 = it->userFloat("PFRelIsoDB04v2");
#endif
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
