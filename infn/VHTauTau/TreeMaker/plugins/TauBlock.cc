#include <iostream>
#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Utilities/General/interface/FileInPath.h"

#include "VHTauTau/TreeMaker/plugins/TauBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

namespace tb 
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref)
  {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}
TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _beamSpotInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _beamSpotCorr(iConfig.getParameter<bool>("beamSpotCorr"))
{
}
TauBlock::~TauBlock() { }
void TauBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneTau = new TClonesArray("vhtm::Tau");
  tree->Branch("Tau", &cloneTau, 32000, 2);
  tree->Branch("nTau", &fnTau, "fnTau/I");
}
void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTau->Clear();
  fnTau = 0;

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_inputTag, taus);
  
  if (taus.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(_vtxInputTag, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    if (_beamSpotCorr) {
      iEvent.getByLabel(_beamSpotInputTag, beamSpot);
    }

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
                                               it != taus->end(); ++it) {
      if (fnTau == kMaxTau) {
	edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << fnTau;
	break;
      }
      tauB = new ((*cloneTau)[fnTau++]) vhtm::Tau();

      // Store Tau variables
      tauB->eta    = it->eta();
      tauB->phi    = it->phi();
      tauB->pt     = it->pt();
      tauB->energy = it->energy();
      tauB->charge = it->charge();
      if (tb::isValidRef(it->leadPFChargedHadrCand()) && tb::isValidRef(it->leadPFChargedHadrCand()->trackRef())) {
	reco::TrackRef trk = it->leadPFChargedHadrCand()->trackRef();
        tauB->leadTrkPt      = trk->pt();
        tauB->leadTrkPtError = trk->ptError();
        tauB->leadTrkEta     = trk->eta();
        tauB->leadTrkPhi     = trk->phi();
        tauB->leadTrkCharge  = trk->charge();

        double trkd0 = trk->d0();
        double trkdz = trk->dz();
        if (_beamSpotCorr) {
          if (beamSpot.isValid()) {
            trkd0 = -(trk->dxy(beamSpot->position()));
            trkdz = trk->dz(beamSpot->position());
          }
          else
            edm::LogError("MuonsBlock") << "Error >> Failed to get BeamSpot for label: "
                                        << _beamSpotInputTag;
        }
        tauB->leadTrkD0      = trkd0;
        tauB->leadTrkD0Error = trk->d0Error();
        tauB->leadTrkDz      = trkdz;
        tauB->leadTrkDzError = trk->dzError();
        if (0) std::cout << trk->pt() << " "
                         << trk->eta() << " "
                         << trk->phi() << " "
                         << trk->charge() << " "
                         << trk->d0() <<  " "
                         << trk->dz() 
                         << std::endl;

        if (primaryVertices.isValid()) {
          edm::LogInfo("TauBlock") << "Total # Primary Vertices: " << primaryVertices->size();

          // IP of leadPFChargedHadrCand wrt event PV
          reco::VertexCollection::const_iterator vit = primaryVertices->begin(); // Highest sumPt vertex
          tauB->dxyPV = trk->dxy(vit->position());
          tauB->dzPV  = trk->dz(vit->position());

          // IP of leadPFChargedHadrCand wrt closest PV
          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDz = 9999.;
          double vertexDxy = 9999.;
          for (reco::VertexCollection::const_iterator vit  = primaryVertices->begin();
                                                      vit != primaryVertices->end(); ++vit) {
            double dxy = trk->dxy(vit->position());
            double dz  = trk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDxy = dxy;
              vertexDz = dz;
            }
          }
          tauB->vtxIndex = indexVtx;
          tauB->vtxDxy   = vertexDxy;
          tauB->vtxDz    = vertexDz;
          if (0) std::cout << indexVtx << " " << vertexDxy << " " << vertexDz << std::endl;
	}
        else {
  	  edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: "
                                    << _vtxInputTag;
        }
      }
      // Leading particle pT
      tauB->leadChargedParticlePt = it->leadPFChargedHadrCand().isNonnull() 
                                      ? it->leadPFChargedHadrCand()->pt(): 0.;
      tauB->leadNeutralParticlePt = it->leadPFNeutralCand().isNonnull() 
                                      ? it->leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = it->leadPFCand().isNonnull() 
                                      ? it->leadPFCand()->et(): 0.;      
      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = it->signalPFCands().size();
      
      tauB->numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = it->isolationPFCands().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = it->isolationPFGammaCandsEtSum();

      // Signal Constituents
      // Charged hadrons
      int indx = 0;
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFChargedHadrCands().begin(); 
    	                                              iCand != it->signalPFChargedHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          tauB->sigChHadCandPt[indx]  = cand.pt();
          tauB->sigChHadCandEta[indx] = cand.eta();
          tauB->sigChHadCandPhi[indx] = cand.phi();
        }
      }  
      // Neutral hadrons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFNeutrHadrCands().begin(); 
    	                                              iCand != it->signalPFNeutrHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->sigNeHadCandPt[indx]  = cand.pt();
          tauB->sigNeHadCandEta[indx] = cand.eta();
          tauB->sigNeHadCandPhi[indx] = cand.phi();
        }
      }  
      // Photons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFGammaCands().begin(); 
    	                                              iCand != it->signalPFGammaCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->sigGammaCandPt[indx]  = cand.pt();
          tauB->sigGammaCandEta[indx] = cand.eta();
          tauB->sigGammaCandPhi[indx] = cand.phi();
        }
      }  
      // Isolation Constituents
      // Charged hadrons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFChargedHadrCands().begin(); 
    	                                              iCand != it->isolationPFChargedHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          tauB->isoChHadCandPt[indx]  = cand.pt();
          tauB->isoChHadCandEta[indx] = cand.eta();
          tauB->isoChHadCandPhi[indx] = cand.phi();
        }
      }  
      
      // Neutral hadrons
      double ptSum = 0;
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFNeutrHadrCands().begin(); 
	                                              iCand != it->isolationPFNeutrHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->isoNeHadCandPt[indx]  = cand.pt();
          tauB->isoNeHadCandEta[indx] = cand.eta();
          tauB->isoNeHadCandPhi[indx] = cand.phi();
        }
	ptSum += std::abs(cand.pt());
      }  
      tauB->ptSumPFNeutralHadronsIsoCone = ptSum;

      // Photons
      indx = 0; // reset
      for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFGammaCands().begin(); 
      	                                              iCand != it->isolationPFGammaCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tauB->isoGammaCandPt[indx]  = cand.pt();
          tauB->isoGammaCandEta[indx] = cand.eta();
          tauB->isoGammaCandPhi[indx] = cand.phi();
        }
      }  

      // tau id. discriminators
      tauB->decayModeFinding = it->tauID("decayModeFinding");

      // discriminators against muons
      tauB->againstMuonLoose      = it->tauID("againstMuonLoose");
      tauB->againstMuonMedium     = it->tauID("againstMuonMedium");
      tauB->againstMuonTight      = it->tauID("againstMuonTight");

      tauB->againstMuonLoose2     = it->tauID("againstMuonLoose2");
      tauB->againstMuonMedium2    = it->tauID("againstMuonMedium2");
      tauB->againstMuonTight2     = it->tauID("againstMuonTight2");

      // discriminators against electrons
      tauB->againstElectronLoose  = it->tauID("againstElectronLoose");
      tauB->againstElectronMedium = it->tauID("againstElectronMedium");
      tauB->againstElectronTight  = it->tauID("againstElectronTight");
      tauB->againstElectronMVA    = it->tauID("againstElectronMVA");

      tauB->againstElectronVLooseMVA2  = it->tauID("againstElectronVLooseMVA2");
      tauB->againstElectronLooseMVA2   = it->tauID("againstElectronLooseMVA2");
      tauB->againstElectronMediumMVA2  = it->tauID("againstElectronMediumMVA2");
      tauB->againstElectronTightMVA2   = it->tauID("againstElectronTightMVA2");
      
      tauB->againstElectronLooseMVA3  = it->tauID("againstElectronLooseMVA3");
      tauB->againstElectronMediumMVA3  = it->tauID("againstElectronMediumMVA3");
      tauB->againstElectronTightMVA3  = it->tauID("againstElectronTightMVA3");
      tauB->againstElectronVTightMVA3  = it->tauID("againstElectronVTightMVA3");

      // Obsolete
      tauB->pfElectronMVA         = it->leadPFCand().isNonnull() 
                                  ? it->leadPFCand()->mva_e_pi() : 1.;

      tauB->byVLooseIsolation  = it->tauID("byVLooseIsolation");
      tauB->byLooseIsolation   = it->tauID("byLooseIsolation");
      tauB->byMediumIsolation  = it->tauID("byMediumIsolation");
      tauB->byTightIsolation   = it->tauID("byTightIsolation");

      // MVA based isolation
      tauB->byLooseIsolationMVA                    = it->tauID("byLooseIsolationMVA"); 
      tauB->byMediumIsolationMVA                   = it->tauID("byMediumIsolationMVA"); 
      tauB->byTightIsolationMVA                    = it->tauID("byTightIsolationMVA"); 

      tauB->byLooseIsolationMVA2                   = it->tauID("byLooseIsolationMVA2"); 
      tauB->byMediumIsolationMVA2                  = it->tauID("byMediumIsolationMVA2"); 
      tauB->byTightIsolationMVA2                   = it->tauID("byTightIsolationMVA2"); 
  
      // DB Corrected Isolation
      tauB->byVLooseCombinedIsolationDeltaBetaCorr = it->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      tauB->byLooseCombinedIsolationDeltaBetaCorr  = it->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      tauB->byMediumCombinedIsolationDeltaBetaCorr = it->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      tauB->byTightCombinedIsolationDeltaBetaCorr  = it->tauID("byTightCombinedIsolationDeltaBetaCorr");

      tauB->byVLooseIsolationDeltaBetaCorr         = it->tauID("byVLooseIsolationDeltaBetaCorr");
      tauB->byLooseIsolationDeltaBetaCorr          = it->tauID("byLooseIsolationDeltaBetaCorr");
      tauB->byMediumIsolationDeltaBetaCorr         = it->tauID("byMediumIsolationDeltaBetaCorr");
      tauB->byTightIsolationDeltaBetaCorr          = it->tauID("byTightIsolationDeltaBetaCorr");

      tauB->byLooseCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      tauB->byMediumCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tauB->byTightCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");

      // kinematic variables for PFJet associated to PFTau
      if (tb::isValidRef(it->pfJetRef())) {
	const reco::PFJetRef &jtr = it->pfJetRef();
        tauB->jetPt  = jtr->pt();
        tauB->jetEta = jtr->eta();
        tauB->jetPhi = jtr->phi();
        if (0) std::cout << jtr->pt() << " " << jtr->eta() << " " << jtr->phi() << std::endl;
      }

      // NEW quantities
      tauB->emFraction              = it->emFraction(); 
      tauB->maximumHCALPFClusterEt  = it->maximumHCALPFClusterEt();
      tauB->ecalStripSumEOverPLead  = it->ecalStripSumEOverPLead();
      tauB->bremsRecoveryEOverPLead = it->bremsRecoveryEOverPLead();
      tauB->hcalTotOverPLead        = it->hcalTotOverPLead();
      tauB->hcalMaxOverPLead        = it->hcalMaxOverPLead();
      tauB->hcal3x3OverPLead        = it->hcal3x3OverPLead();

      tauB->etaetaMoment = it->etaetaMoment();
      tauB->phiphiMoment = it->phiphiMoment();
      tauB->phiphiMoment = it->etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      tauB->vx = vertex.x();             
      tauB->vy = vertex.y();             
      tauB->vz = vertex.z();             

      tauB->zvertex = it->vz(); // distance from the primary vertex
      tauB->mass    = it->p4().M();
      tauB->ltsipt  = TMath::Abs(it->leadPFChargedHadrCandsignedSipt());
    }
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: " 
                              << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
