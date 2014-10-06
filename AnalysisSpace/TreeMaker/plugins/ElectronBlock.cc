#include <iostream> 
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVector3.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/ElectronBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", false)),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  electronTag_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc", edm::InputTag("selectedPatElectrons"))),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  electronToken_(consumes<pat::ElectronCollection>(electronTag_))
{
}
ElectronBlock::~ElectronBlock() {
}
void ElectronBlock::beginJob()
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::Electron>();
  tree->Branch("Electron", "std::vector<vhtm::Electron>", &list_, 32000, 2);
  tree->Branch("nElectron", &fnElectron_, "fnElectron_/I");
}
void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnElectron_ = 0;

  edm::Handle<pat::ElectronCollection> electrons;
  bool found = iEvent.getByToken(electronToken_, electrons);

  if (found && electrons.isValid()) {
    edm::Handle<reco::BeamSpot> beamSpot;
    if (bsCorr_) iEvent.getByToken(bsToken_, beamSpot);

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (auto it = electrons->begin(); it != electrons->end(); ++it) {
      if (fnElectron_ == kMaxElectron_) {
        edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = "
                                      << fnElectron_;
        break;
      }
      bool hasGsfTrack = it->gsfTrack().isNonnull() ? true : false;

      vhtm::Electron electron;
      electron.ecalDriven      = it->ecalDrivenSeed();
      electron.eta             = it->eta();
      electron.phi             = it->phi();
      electron.pt              = it->pt();
      electron.hasGsfTrack     = hasGsfTrack;
      electron.energy          = it->energy();
      electron.caloEnergy      = it->ecalEnergy();
      electron.caloEnergyError = it->ecalEnergyError();
      electron.charge          = it->charge();
  
      if (hasGsfTrack) {
        reco::GsfTrackRef tk = it->gsfTrack();
        electron.trackPt      = tk->pt();
        electron.trackPtError = tk->ptError();

        // Hit pattern
        const reco::HitPattern& hitp = tk->hitPattern();
        electron.pixHits = hitp.numberOfValidPixelHits();
        electron.trkHits = hitp.numberOfValidTrackerHits();

        electron.nValidHits  = tk->numberOfValidHits();
        electron.missingHits = tk->trackerExpectedHitsInner().numberOfHits();

        double trkd0 = tk->d0();
        double trkdz = tk->dz();
        if (bsCorr_) {
          if (beamSpot.isValid()) {
            trkd0 = -(tk->dxy(beamSpot->position()));
            trkdz = tk->dz(beamSpot->position());
          }
          else
            edm::LogError("ElectronBlock") << "Error >> Failed to get BeamSpot for label: "
                                           << bsTag_;
        }
        electron.trkD0      = trkd0;
        electron.trkD0Error = tk->d0Error();
        electron.trkDz      = trkdz;
        electron.trkDzError = tk->dzError();

        if (primaryVertices.isValid()) {
          auto vit = primaryVertices->begin(); // Highest sumPt vertex
          electron.dxyPV = tk->dxy(vit->position());
          electron.dzPV = tk->dz(vit->position());

          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDistZ = 9999.;
          edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
          for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
            double dxy = tk->dxy(vit->position());
            double dz  = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy, 2) + pow(dz, 2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
          electron.vtxDist3D = minVtxDist3D;
          electron.vtxIndex = indexVtx;
          electron.vtxDistZ = vertexDistZ;
        }
        else {
          edm::LogError("ElectronBlock") << "Error >> Failed to get VertexCollection for label: "
                                         << vertexTag_;
        }
      }
      // ID variables
      electron.hoe           = it->hcalOverEcal();
      electron.hoeDepth1     = it->hcalDepth1OverEcal();
      electron.eop           = it->eSuperClusterOverP();
      electron.sigmaEtaEta   = it->sigmaEtaEta();
      electron.sigmaIEtaIEta = it->sigmaIetaIeta();
      electron.deltaPhiTrkSC = it->deltaPhiSuperClusterTrackAtVtx();
      electron.deltaEtaTrkSC = it->deltaEtaSuperClusterTrackAtVtx();
      electron.classif       = it->classification();
      electron.e1x5overe5x5  = (it->e5x5() > 0) ? it->e1x5()/it->e5x5() : 0;
      electron.e2x5overe5x5  = (it->e5x5() > 0) ? it->e2x5Max()/it->e5x5() : 0;

      // Iso variables
      electron.isoEcal03     = it->dr03EcalRecHitSumEt();
      electron.isoHcal03     = it->dr03HcalTowerSumEt();
      electron.isoTrk03      = it->dr03TkSumPt();
      electron.isoEcal04     = it->dr04EcalRecHitSumEt(); // ecalIso
      electron.isoHcal04     = it->dr04HcalTowerSumEt();  // hcalIso
      electron.isoTrk04      = it->dr04TkSumPt();         // trackIso
      electron.isoRel03      = (it->dr03EcalRecHitSumEt()
                              + it->dr03HcalTowerSumEt()
                              + it->dr03TkSumPt())/it->pt();
      electron.isoRel04      = (it->dr04EcalRecHitSumEt()
                              + it->dr04HcalTowerSumEt()
                              + it->dr04TkSumPt())/it->pt();

      // SC associated with electron
      electron.scEn  = it->superCluster()->energy();
      electron.scEta = it->superCluster()->eta();
      electron.scPhi = it->superCluster()->phi();
      electron.scET  = it->superCluster()->energy()/cosh(it->superCluster()->eta());
      electron.scRawEnergy = it->superCluster()->rawEnergy();

      electron.dist_vec       = it->userFloat("dist");
      electron.dCotTheta      = it->userFloat("dcot");
      electron.hasMatchedConv = (it->userInt("antiConv") ? false : true);

      electron.relIso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      // PF Isolation
      reco::GsfElectron::PflowIsolationVariables pfIso = it->pfIsolationVariables();
      electron.sumChargedHadronPt = pfIso.sumChargedHadronPt;
      electron.sumPUPt = pfIso.sumPUPt;
      float absiso = pfIso.sumChargedHadronPt + std::max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt);
      float iso = absiso/(it->p4().pt());
      electron.pfRelIso = iso;

#if 0
      // PF based isolation
      electron.pfRelIso = it->userFloat("PFRelIsoDB04");
#endif
      // PFlow isolation information
      electron.chargedHadronIso = it->chargedHadronIso();
      electron.neutralHadronIso = it->neutralHadronIso();
      electron.photonIso        = it->photonIso();
  
      // IP information
      electron.dB  = it->dB(pat::Electron::PV2D);
      electron.edB = it->edB(pat::Electron::PV2D);

      electron.dB3d  = it->dB(pat::Electron::PV3D);
      electron.edB3d = it->edB(pat::Electron::PV3D);

      // Bremstrahlung information
      electron.nBrems = it->numberOfBrems();
      electron.fbrem  = it->fbrem();

#if 0
      // MVA
      electron.mva               = it->userFloat("mva");
      electron.mvaPOGTrig        = it->userFloat("mvaPOGTrig");
      electron.mvaPOGNonTrig     = it->userFloat("mvaPOGNonTrig");
#endif
      electron.mvaPreselection   = (it->userInt("mvaPreselection") ? true : false);
      electron.isTriggerElectron = (it->userInt("isTriggerElectron") ? true : false);
#if 0
      // MVA Iso
      electron.isoMVA = it->userFloat("eleIsoMVA");
#endif
      // Fiducial flag
      int fidFlag = 0;
      if (it->isEB()) fidFlag |= (1 << 0);
      if (it->isEE()) fidFlag |= (1 << 1);
      if (it->isEBEtaGap()) fidFlag |= (1 << 2);
      if (it->isEBPhiGap()) fidFlag |= (1 << 3);
      if (it->isEERingGap()) fidFlag |= (1 << 4);
      if (it->isEEDeeGap()) fidFlag |= (1 << 5);
      if (it->isEBEEGap()) fidFlag |= (1 << 6);
      electron.fidFlag = fidFlag;

      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      electron.vx = vertex.x();
      electron.vy = vertex.y();
      electron.vz = vertex.z();
     
#if 0
      electron.pfRelIso03v1   = it->userFloat("PFRelIso03v1");
      electron.pfRelIso03v2   = it->userFloat("PFRelIso03v2");
      electron.pfRelIsoDB03v1 = it->userFloat("PFRelIsoDB03v1");
      electron.pfRelIsoDB03v2 = it->userFloat("PFRelIsoDB03v2");
      electron.pfRelIsoDB03v3 = it->userFloat("PFRelIsoDB03v3");

      electron.pfRelIso04v1   = it->userFloat("PFRelIso04v1");
      electron.pfRelIso04v2   = it->userFloat("PFRelIso04v2");
      electron.pfRelIsoDB04v1 = it->userFloat("PFRelIsoDB04v1");
      electron.pfRelIsoDB04v2 = it->userFloat("PFRelIsoDB04v2");
      electron.pfRelIsoDB04v3 = it->userFloat("PFRelIsoDB04v3");

      // 2012
      electron.pfRelIso03   = it->userFloat("PFRelIso03");
      electron.pfRelIso04   = it->userFloat("PFRelIso04");
      electron.pfRelIsoDB03 = it->userFloat("PFRelIsoDB03");
      electron.pfRelIsoDB04 = it->userFloat("PFRelIsoDB04");
#endif
      list_->push_back(electron);
    }
    fnElectron_ = list_->size();
  }
  else {
    edm::LogError("ElectronBlock") << "Error >> Failed to get pat::Electron Collection for label: "
                                   << electronTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronBlock);
