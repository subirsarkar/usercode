#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "VHTauTau/TreeMaker/plugins/ElectronBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

// Constructor
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _bsInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc"))
{
}
ElectronBlock::~ElectronBlock() { 
}
void ElectronBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneElectron = new TClonesArray("vhtm::Electron");
  tree->Branch("Electron", &cloneElectron, 32000, 2);
  tree->Branch("nElectron", &fnElectron,  "fnElectron/I");
}
void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneElectron->Clear();
  fnElectron = 0;

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(_bsInputTag, beamSpot);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_electronInputTag, electrons);

  if (electrons.isValid()) {
    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (std::vector<pat::Electron>::const_iterator it  = electrons->begin(); 
                                                    it != electrons->end(); ++it) {
      if (fnElectron == kMaxElectron) {
	edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = " 
                                      << fnElectron; 
	break;
      }
      bool hasGsfTrack  = it->gsfTrack().isNonnull() ? true : false;
      reco::GsfTrackRef tk = it->gsfTrack();

      electronB = new ((*cloneElectron)[fnElectron++]) vhtm::Electron();
      electronB->ecalDriven = it->ecalDrivenSeed();
      electronB->eta         = it->eta();
      electronB->phi         = it->phi();
      electronB->pt          = it->pt();
      electronB->hasGsfTrack = hasGsfTrack;
      electronB->energy      = it->energy();
      electronB->caloEnergy  = it->ecalEnergy();
      electronB->caloEnergyError = it->ecalEnergyError();
      electronB->charge      = it->charge();
  
      if (hasGsfTrack) {
        electronB->trackPt      = tk->pt();
        electronB->trackPtError = tk->ptError();

	// Hit pattern
	const reco::HitPattern& hitp = tk->hitPattern();
	electronB->pixHits = hitp.numberOfValidPixelHits();
	electronB->trkHits = hitp.numberOfValidTrackerHits();

        electronB->nValidHits   = tk->numberOfValidHits(); 
        electronB->missingHits  = tk->trackerExpectedHitsInner().numberOfHits();

        electronB->trkD0        = tk->d0();
        electronB->trkD0Error   = tk->d0Error();
      }
      // ID variables
      electronB->hoe           = it->hcalOverEcal();
      electronB->hoeDepth1     = it->hcalDepth1OverEcal();
      electronB->eop           = it->eSuperClusterOverP(); 
      electronB->sigmaEtaEta   = it->sigmaEtaEta();
      electronB->sigmaIEtaIEta = it->sigmaIetaIeta();
      electronB->deltaPhiTrkSC = it->deltaPhiSuperClusterTrackAtVtx();
      electronB->deltaEtaTrkSC = it->deltaEtaSuperClusterTrackAtVtx();
      electronB->classif       = it->classification();
      electronB->e1x5overe5x5  = (it->e5x5() > 0) ? (it->e1x5()/it->e5x5()) : 0;
      electronB->e2x5overe5x5  = (it->e5x5() > 0) ? (it->e2x5Max()/it->e5x5()) : 0;

      // Iso variables
      electronB->isoEcal03 = it->dr03EcalRecHitSumEt();
      electronB->isoHcal03 = it->dr03HcalTowerSumEt();
      electronB->isoTrk03  = it->dr03TkSumPt();
      electronB->isoEcal04 = it->dr04EcalRecHitSumEt(); // ecalIso
      electronB->isoHcal04 = it->dr04HcalTowerSumEt(); // hcalIso
      electronB->isoTrk04  = it->dr04TkSumPt(); // trackIso
      electronB->isoRel03  = (it->dr03EcalRecHitSumEt()
                            + it->dr03HcalTowerSumEt()
                            + it->dr03TkSumPt())/it->pt();
      electronB->isoRel04  = (it->dr04EcalRecHitSumEt()
                            + it->dr04HcalTowerSumEt()
                            + it->dr04TkSumPt())/it->pt();

      // SC associated with electron
      electronB->scEn        = it->superCluster()->energy();
      electronB->scEta       = it->superCluster()->eta();
      electronB->scPhi       = it->superCluster()->phi();
      electronB->scET        = it->superCluster()->energy()/cosh(it->superCluster()->eta());
      electronB->scRawEnergy = it->superCluster()->rawEnergy();

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (hasGsfTrack) {
        if (primaryVertices.isValid()) {
	  edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
          for (reco::VertexCollection::const_iterator vit  = primaryVertices->begin(); 
                                                      vit != primaryVertices->end(); ++vit) {
            double dxy = tk->dxy(vit->position());
            double dz  = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy, 2) + pow(dz, 2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
        } 
        else {
	  edm::LogError("ElectronBlock") << "Error >> Failed to get VertexCollection for label: " 
                                         << _vtxInputTag;
        }      
      }
      // Vertex association variables
      electronB->vtxDist3D = minVtxDist3D;
      electronB->vtxIndex  = indexVtx;
      electronB->vtxDistZ  = vertexDistZ;

      electronB->relIso   = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      // PF based isolation
      electronB->pfRelIso = 
	it->chargedHadronIso() + std::max(0.0, it->neutralHadronIso() + it->photonIso() - 0.5 * it->puChargedHadronIso())/it->pt();

      // PFlow isolation information
      electronB->chargedHadronIso = it->chargedHadronIso();
      electronB->neutralHadronIso = it->neutralHadronIso();
      electronB->photonIso        = it->photonIso();
  
      // IP information
      electronB->dB  = it->dB(pat::Electron::PV2D);
      electronB->edB = it->edB(pat::Electron::PV2D);

      electronB->dB3d  = it->dB(pat::Electron::PV3D);
      electronB->edB3d = it->edB(pat::Electron::PV3D);

      // Bremstrahlung information
      electronB->nBrems = it->numberOfBrems();
      electronB->fbrem  = it->fbrem();

      electronB->idMVA = it->electronID("mvaNonTrigV0");
 
      // Fiducial flag 
      int fidFlag = 0;
      if (it->isEB())        fidFlag |= (1 << 0);
      if (it->isEE())        fidFlag |= (1 << 1);
      if (it->isEBEtaGap())  fidFlag |= (1 << 2);
      if (it->isEBPhiGap())  fidFlag |= (1 << 3);
      if (it->isEERingGap()) fidFlag |= (1 << 4);
      if (it->isEEDeeGap())  fidFlag |= (1 << 5);
      if (it->isEBEEGap())   fidFlag |= (1 << 6);
      electronB->fidFlag = fidFlag;

      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      electronB->vx = vertex.x();             
      electronB->vy = vertex.y();             
      electronB->vz = vertex.z();             
    }
  } 
  else {
    edm::LogError("ElectronBlock") << "Error >> Failed to get pat::Electron Collection for label: " 
                                   << _electronInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronBlock);
