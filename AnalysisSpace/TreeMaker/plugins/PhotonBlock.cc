#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "AnalysisSpace/TreeMaker/plugins/PhotonBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
PhotonBlock::PhotonBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  photonTag_(iConfig.getParameter<edm::InputTag>("photonSrc")),
  photonToken_(consumes<pat::PhotonCollection>(photonTag_))
{}
PhotonBlock::~PhotonBlock() {
  delete list_;
}
void PhotonBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::Photon>();
  tree->Branch("Photon", "std::vector<vhtm::Photon>", &list_, 32000, 2);
  tree->Branch("nPhoton", &fnPhoton_, "fnPhoton_/I");
}
void PhotonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  list_->clear();
  fnPhoton_ = 0;

  edm::Handle<pat::PhotonCollection> photons;
  bool found = iEvent.getByToken(photonToken_, photons);

  if (found && photons.isValid()) {
    edm::LogInfo("PhotonBlock") << "Total # PAT Photons: " << photons->size();
    for (pat::Photon const& v: *photons) {
      if (list_->size() == kMaxPhoton) {
	edm::LogInfo("PhotonBlock") << "Too many PAT Photon, fnPhoton = " 
                                    << fnPhoton_; 
	break;
      }
      vhtm::Photon photon;
      photon.et     = v.et();
      photon.eta    = v.eta();
      photon.clusterEta = v.caloPosition().eta();
      photon.phi    = v.phi();
      photon.clusterPhi = v.caloPosition().phi();
      photon.energy = v.energy();
      photon.theta  = v.theta();
      photon.vx     = v.vx();
      photon.vy     = v.vy();
      photon.vz     = v.vz();

      const reco::SuperClusterRef sCluster = v.superCluster(); 
      if (sCluster.isNonnull()) {
	photon.scEnergy    = sCluster->energy();
	photon.scEta       = sCluster->eta();
	photon.scPhi       = sCluster->phi();
	photon.scSize      = sCluster->clustersSize();
	photon.scEtaWidth  = sCluster->etaWidth();
	photon.scPhiWidth  = sCluster->phiWidth();
	photon.scEt        = sCluster->energy()/cosh(sCluster->eta());
	photon.scRawEnergy = sCluster->rawEnergy();
	photon.scx         = sCluster->x();
	photon.scy         = sCluster->y();
	photon.scz         = sCluster->z();
      }
      photon.isoEcalRecHit03    = v.ecalRecHitSumEtConeDR03();
      photon.isoHcalRecHit03    = v.hcalTowerSumEtConeDR03();
      photon.isoSolidTrkCone03  = v.trkSumPtSolidConeDR03();
      photon.isoHollowTrkCone03 = v.trkSumPtHollowConeDR03();
      photon.nTrkSolidCone03    = v.nTrkSolidConeDR03();
      photon.nTrkHollowCone03   = v.nTrkHollowConeDR03();

      photon.isoEcalRecHit04    = v.ecalRecHitSumEtConeDR04();
      photon.isoHcalRecHit04    = v.hcalTowerSumEtConeDR04();
      photon.isoSolidTrkCone04  = v.trkSumPtSolidConeDR04();
      photon.isoHollowTrkCone04 = v.trkSumPtHollowConeDR04();
      photon.nTrkSolidCone04    = v.nTrkSolidConeDR04();
      photon.nTrkHollowCone04   = v.nTrkHollowConeDR04();

      photon.hasPixelSeed       = v.hasPixelSeed(); 
      photon.passElectronVeto   = v.passElectronVeto();

      photon.chargedHadIso      = v.chargedHadronIso();
      photon.neutralHadIso      = v.neutralHadronIso();
      photon.photonIso          = v.photonIso();
      photon.puChargedHadIso    = v.puChargedHadronIso();

      int fidFlag = 0;
      if (v.isEB())        fidFlag |= (1 << 0);
      if (v.isEE())        fidFlag |= (1 << 1);
      if (v.isEBEtaGap())  fidFlag |= (1 << 2);
      if (v.isEBPhiGap())  fidFlag |= (1 << 3);
      if (v.isEERingGap()) fidFlag |= (1 << 4);
      if (v.isEEDeeGap())  fidFlag |= (1 << 5);
      if (v.isEBEEGap())   fidFlag |= (1 << 6);
      photon.fidFlag = fidFlag;

      photon.isEB               = v.isEB() ? true : false;
      photon.isEE               = v.isEE() ? true : false;
      photon.isEBGap            = v.isEBGap() ? true : false ;
      photon.isEEGap            = v.isEEGap() ? true : false;
      photon.isEBEEGap          = v.isEBEEGap() ? true : false;

      photon.r9                 = v.r9();
      photon.hoe                = v.hadronicOverEm();
      photon.sigmaEtaEta        = v.sigmaEtaEta();
      photon.sigmaIEtaIEta      = v.sigmaIetaIeta();
      photon.e1x5               = v.e1x5();
      photon.e2x5               = v.e2x5();
      photon.e3x3               = v.e3x3();
      photon.e5x5               = v.e5x5();
      photon.r1x5               = v.r1x5();
      photon.r2x5               = v.r2x5();
      photon.maxEnergyXtal      = v.maxEnergyXtal();

      for (const pat::Photon::IdPair& pa: v.photonIDs())
	photon.idmap[pa.first] = pa.second;

      photon.hasConversionTracks = v.hasConversionTracks();      
      if (v.hasConversionTracks()) {
        const reco::ConversionRefVector conversions = v.conversions();
        for (edm::RefVector<reco::ConversionCollection>::const_iterator jt  = conversions.begin();
                                                                        jt != conversions.end(); 
                                                                      ++jt) 
	{
          const reco::Conversion& obj = (**jt);
	  if (obj.nTracks() < 2 or
              !obj.conversionVertex().isValid()) continue;
          photon.nTracks = obj.nTracks();
          photon.isConverted = obj.isConverted();
          photon.pairInvMass = obj.pairInvariantMass();
          photon.pairCotThetaSeparation = obj.pairCotThetaSeparation();

	  math::XYZVectorF  mom = obj.pairMomentum();
          photon.pairPx = mom.x();
          photon.pairPy = mom.y();
          photon.pairPz = mom.z();

          const reco::Vertex &vertex = obj.conversionVertex();
          photon.conv_vx = vertex.x();
          photon.conv_vy = vertex.y();
          photon.conv_vz = vertex.z();

	  photon.eovp              = obj.EoverP();
	  photon.zpv               = obj.zOfPrimaryVertexFromTracks();
	  photon.distOfMinApproach = obj.distOfMinimumApproach();
	  photon.dPhiTracksAtVtx   = obj.dPhiTracksAtVtx();
	  photon.dPhiTracksAtEcal  = obj.dPhiTracksAtEcal();
	  photon.dEtaTracksAtEcal  = obj.dEtaTracksAtEcal();
        }    
      }
      list_->push_back(photon);
    }
    fnPhoton_ = list_->size();
  }
  else {
    edm::LogError("PhotonBlock") << "Error >> Failed to get pat::Photon for label: " 
                                 << photonTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonBlock);

