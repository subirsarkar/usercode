#include <iostream>

#include "TTree.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "Utilities/General/interface/FileInPath.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/TauBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

namespace tb
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref) {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}
TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  tauTag_(iConfig.getUntrackedParameter<edm::InputTag>("patTauSrc", edm::InputTag("selectedPatTaus"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  tauToken_(consumes<pat::TauCollection>(tauTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_))
{
}
TauBlock::~TauBlock() { }
void TauBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Tau>();
  tree->Branch("Tau", "std::vector<vhtm::Tau>", &list_, 32000, 2);
  tree->Branch("nTau", &fnTau_, "fnTau_/I");
}
void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTau_ = 0;

  edm::Handle<pat::TauCollection> taus;
  bool found = iEvent.getByToken(tauToken_, taus);
  
  if (found && taus.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    if (bsCorr_)
      iEvent.getByToken(bsToken_, beamSpot);

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    for (auto it = taus->begin(); it != taus->end(); ++it) {
      if (fnTau_ == kMaxTau_) {
        edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << fnTau_;
        break;
      }
      vhtm::Tau tau;

      // Store Tau variables
      tau.eta    = it->eta();
      tau.phi    = it->phi();
      tau.pt     = it->pt();
      tau.energy = it->energy();
      tau.charge = it->charge();
      if (it->leadPFChargedHadrCand().isNonnull() && tb::isValidRef(it->leadPFChargedHadrCand()->trackRef())) {
        reco::TrackRef trk = it->leadPFChargedHadrCand()->trackRef();
        tau.leadTrkPt      = trk->pt();
        tau.leadTrkPtError = trk->ptError();
        tau.leadTrkEta     = trk->eta();
        tau.leadTrkPhi     = trk->phi();
        tau.leadTrkCharge  = trk->charge();

        double trkd0 = trk->d0();
        double trkdz = trk->dz();
        if (bsCorr_) {
          if (beamSpot.isValid()) {
            trkd0 = -(trk->dxy(beamSpot->position()));
            trkdz = trk->dz(beamSpot->position());
          }
          else
            edm::LogError("TauBlock") << "Error >> Failed to get reco::BeamSpot for label: "
                                        << bsTag_;
        }
        tau.leadTrkD0      = trkd0;
        tau.leadTrkD0Error = trk->d0Error();
        tau.leadTrkDz      = trkdz;
        tau.leadTrkDzError = trk->dzError();

        if (primaryVertices.isValid()) {
          edm::LogInfo("TauBlock") << "Total # Primary Vertices: " 
                                   << primaryVertices->size();

          // IP of leadPFChargedHadrCand wrt event PV
          auto vit = primaryVertices->begin(); // Highest sumPt vertex
          tau.dxyPV = trk->dxy(vit->position());
          tau.dzPV  = trk->dz(vit->position());

          // IP of leadPFChargedHadrCand wrt closest PV
          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDz = 9999.;
          double vertexDxy = 9999.;
          for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
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
          tau.vtxIndex = indexVtx;
          tau.vtxDxy   = vertexDxy;
          tau.vtxDz    = vertexDz;
        }
        else {
          edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: "
                                    << vertexTag_;
        }
      }
      // Leading particle pT
      tau.leadChargedParticlePt = 
           it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->pt(): 0.;
      tau.leadNeutralParticlePt = 
           it->leadPFNeutralCand().isNonnull() ? it->leadPFNeutralCand()->et(): 0.;
      tau.leadParticlePt = 
           it->leadPFCand().isNonnull() ? it->leadPFCand()->et(): 0.;

      // Number of charged/neutral candidates and photons in different cones
      tau.numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tau.numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tau.numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tau.numParticlesSignalCone      = it->signalPFCands().size();
      
      tau.numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tau.numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tau.numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tau.numParticlesIsoCone      = it->isolationPFCands().size();
      
      // Signal Constituents
      // Charged hadrons
      int indx = 0;
      for (auto iCand  = it->signalPFChargedHadrCands().begin();
                iCand != it->signalPFChargedHadrCands().end(); ++iCand,indx++) {
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          const reco::PFCandidate& cand = (**iCand);
          tau.sigChHadCandPt[indx]  = cand.pt();
          tau.sigChHadCandEta[indx] = cand.eta();
          tau.sigChHadCandPhi[indx] = cand.phi();
        }
      }
      // Neutral hadrons
      indx = 0; // reset
      for (auto iCand  = it->signalPFNeutrHadrCands().begin();
                iCand != it->signalPFNeutrHadrCands().end(); ++iCand,indx++) {
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          const reco::PFCandidate& cand = (**iCand);
          tau.sigNeHadCandPt[indx]  = cand.pt();
          tau.sigNeHadCandEta[indx] = cand.eta();
          tau.sigNeHadCandPhi[indx] = cand.phi();
        }
      }
      // Photons
      indx = 0; // reset
      for (auto iCand  = it->signalPFGammaCands().begin();
                iCand != it->signalPFGammaCands().end(); ++iCand,indx++) {
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          const reco::PFCandidate& cand = (**iCand);
          tau.sigGammaCandPt[indx]  = cand.pt();
          tau.sigGammaCandEta[indx] = cand.eta();
          tau.sigGammaCandPhi[indx] = cand.phi();
        }
      }
      // Isolation Constituents
      // Charged hadrons
      indx = 0; // reset
      for (auto iCand  = it->isolationPFChargedHadrCands().begin();
                iCand != it->isolationPFChargedHadrCands().end(); ++iCand,indx++) {
        if (indx < vhtm::Tau::kMaxPFChargedCand) {
          const reco::PFCandidate& cand = (**iCand);
          tau.isoChHadCandPt[indx]  = cand.pt();
          tau.isoChHadCandEta[indx] = cand.eta();
          tau.isoChHadCandPhi[indx] = cand.phi();
        }
      }
      
      // Neutral hadrons
      double ptSum = 0;
      indx = 0; // reset
      for (auto iCand  = it->isolationPFNeutrHadrCands().begin();
                iCand != it->isolationPFNeutrHadrCands().end(); ++iCand,indx++) {
        const reco::PFCandidate& cand = (**iCand);
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          tau.isoNeHadCandPt[indx]  = cand.pt();
          tau.isoNeHadCandEta[indx] = cand.eta();
          tau.isoNeHadCandPhi[indx] = cand.phi();
        }
        ptSum += std::abs(cand.pt());
      }
      // Photons
      indx = 0; // reset
      for (auto iCand  = it->isolationPFGammaCands().begin();
                iCand != it->isolationPFGammaCands().end(); ++iCand,indx++) {
        if (indx < vhtm::Tau::kMaxPFNeutralCand) {
          const reco::PFCandidate& cand = (**iCand);
          tau.isoGammaCandPt[indx]  = cand.pt();
          tau.isoGammaCandEta[indx] = cand.eta();
          tau.isoGammaCandPhi[indx] = cand.phi();
        }
      }
      // PtSum
      tau.ptSumChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tau.ptSumPhotonsIsoCone        = it->isolationPFGammaCandsEtSum();
      tau.ptSumNeutralHadronsIsoCone = ptSum;

      // tau id. discriminators
      tau.decayModeFinding = it->tauID("decayModeFinding");
      tau.decayModeFindingNewDMs = it->tauID("decayModeFindingNewDMs");
      tau.decayModeFindingOldDMs = it->tauID("decayModeFindingOldDMs");

      // discriminators against muons
      tau.againstMuonLoose  = it->tauID("againstMuonLoose");
      tau.againstMuonMedium = it->tauID("againstMuonMedium");
      tau.againstMuonTight  = it->tauID("againstMuonTight");

      tau.againstMuonLoose2  = it->tauID("againstMuonLoose2");
      tau.againstMuonMedium2 = it->tauID("againstMuonMedium2");
      tau.againstMuonTight2  = it->tauID("againstMuonTight2");

      tau.againstMuonLoose3  = it->tauID("againstMuonLoose3");
      tau.againstMuonTight3 = it->tauID("againstMuonTight3");

      // discriminators against electrons
      tau.againstElectronLoose  = it->tauID("againstElectronLoose");
      tau.againstElectronMedium = it->tauID("againstElectronMedium");
      tau.againstElectronTight  = it->tauID("againstElectronTight");

      tau.againstElectronLooseMVA5  = it->tauID("againstElectronLooseMVA5");
      tau.againstElectronMediumMVA5 = it->tauID("againstElectronMediumMVA5");
      tau.againstElectronTightMVA5  = it->tauID("againstElectronTightMVA5");

#if 0
      tau.againstElectronMVA    = it->tauID("againstElectronMVA");

      tau.againstElectronVLooseMVA2 = it->tauID("againstElectronVLooseMVA2");
      tau.againstElectronLooseMVA2  = it->tauID("againstElectronLooseMVA2");
      tau.againstElectronMediumMVA2 = it->tauID("againstElectronMediumMVA2");
      tau.againstElectronTightMVA2  = it->tauID("againstElectronTightMVA2");
      
      tau.againstElectronLooseMVA3  = it->tauID("againstElectronLooseMVA3");
      tau.againstElectronMediumMVA3 = it->tauID("againstElectronMediumMVA3");
      tau.againstElectronTightMVA3  = it->tauID("againstElectronTightMVA3");
      tau.againstElectronVTightMVA3 = it->tauID("againstElectronVTightMVA3");

      // Obsolete
      tau.pfElectronMVA = it->leadPFCand().isNonnull()
                        ? it->leadPFCand()->mva_e_pi() : 1.;
      tau.byVLooseIsolation = it->tauID("byVLooseIsolation");
      tau.byLooseIsolation  = it->tauID("byLooseIsolation");
      tau.byMediumIsolation = it->tauID("byMediumIsolation");
      tau.byTightIsolation  = it->tauID("byTightIsolation");

      // MVA based isolation
      tau.byLooseIsolationMVA  = it->tauID("byLooseIsolationMVA");
      tau.byMediumIsolationMVA = it->tauID("byMediumIsolationMVA");
      tau.byTightIsolationMVA  = it->tauID("byTightIsolationMVA");

      tau.byLooseIsolationMVA2  = it->tauID("byLooseIsolationMVA2");
      tau.byMediumIsolationMVA2 = it->tauID("byMediumIsolationMVA2");
      tau.byTightIsolationMVA2  = it->tauID("byTightIsolationMVA2");
  
#endif
      // DB Corrected Isolation
      tau.byVLooseCombinedIsolationDeltaBetaCorr = it->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      tau.byLooseCombinedIsolationDeltaBetaCorr  = it->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      tau.byMediumCombinedIsolationDeltaBetaCorr = it->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      tau.byTightCombinedIsolationDeltaBetaCorr  = it->tauID("byTightCombinedIsolationDeltaBetaCorr");

#if 0
      tau.byVLooseIsolationDeltaBetaCorr = it->tauID("byVLooseIsolationDeltaBetaCorr");
      tau.byLooseIsolationDeltaBetaCorr  = it->tauID("byLooseIsolationDeltaBetaCorr");
      tau.byMediumIsolationDeltaBetaCorr = it->tauID("byMediumIsolationDeltaBetaCorr");
      tau.byTightIsolationDeltaBetaCorr  = it->tauID("byTightIsolationDeltaBetaCorr");
#endif
      tau.byLooseCombinedIsolationDeltaBetaCorr3Hits  = it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      tau.byMediumCombinedIsolationDeltaBetaCorr3Hits = it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tau.byTightCombinedIsolationDeltaBetaCorr3Hits  = it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      tau.byCombinedIsolationDeltaBetaCorrRaw3Hits    = it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

#if 1
      // Isolation variables suggested by Christian
      tau.chargedIsoPtSum = it->tauID("chargedIsoPtSum");
      tau.neutralIsoPtSum = it->tauID("neutralIsoPtSum");
      tau.puCorrPtSum     = it->tauID("puCorrPtSum");
#endif
#if 0
The available IDs are: 'againstElectronDeadECAL' 'againstElectronLoose' 'againstElectronLooseMVA5' 'againstElectronMVA5category' 
'againstElectronMVA5raw' 'againstElectronMedium' 'againstElectronMediumMVA5' 'againstElectronTight' 'againstElectronTightMVA5' 
'againstElectronVLooseMVA5' 'againstElectronVTightMVA5' 'againstMuonLoose' 'againstMuonLoose2' 'againstMuonLoose3' 'againstMuonLooseMVA' 
'againstMuonMVAraw' 'againstMuonMedium' 'againstMuonMedium2' 'againstMuonMediumMVA' 'againstMuonTight' 'againstMuonTight2' 'againstMuonTight3' 
'againstMuonTightMVA' 'byCombinedIsolationDeltaBetaCorrRaw' 'byCombinedIsolationDeltaBetaCorrRaw3Hits' 'byIsolationMVA3newDMwLTraw' 
'byIsolationMVA3newDMwoLTraw' 'byIsolationMVA3oldDMwLTraw' 'byIsolationMVA3oldDMwoLTraw' 'byLooseCombinedIsolationDeltaBetaCorr' 
'byLooseCombinedIsolationDeltaBetaCorr3Hits' 'byLooseIsolation' 'byLooseIsolationMVA3newDMwLT' 'byLooseIsolationMVA3newDMwoLT' 
'byLooseIsolationMVA3oldDMwLT' 'byLooseIsolationMVA3oldDMwoLT' 'byMediumCombinedIsolationDeltaBetaCorr' 
'byMediumCombinedIsolationDeltaBetaCorr3Hits' 'byMediumIsolationMVA3newDMwLT' 'byMediumIsolationMVA3newDMwoLT' 'byMediumIsolationMVA3oldDMwLT' 
'byMediumIsolationMVA3oldDMwoLT' 'byTightCombinedIsolationDeltaBetaCorr' 'byTightCombinedIsolationDeltaBetaCorr3Hits' 
'byTightIsolationMVA3newDMwLT' 'byTightIsolationMVA3newDMwoLT' 'byTightIsolationMVA3oldDMwLT' 'byTightIsolationMVA3oldDMwoLT' 
'byVLooseCombinedIsolationDeltaBetaCorr' 'byVLooseIsolationMVA3newDMwLT' 'byVLooseIsolationMVA3newDMwoLT' 'byVLooseIsolationMVA3oldDMwLT' 
'byVLooseIsolationMVA3oldDMwoLT' 'byVTightIsolationMVA3newDMwLT' 'byVTightIsolationMVA3newDMwoLT' 'byVTightIsolationMVA3oldDMwLT' 
'byVTightIsolationMVA3oldDMwoLT' 'byVVTightIsolationMVA3newDMwLT' 'byVVTightIsolationMVA3newDMwoLT' 'byVVTightIsolationMVA3oldDMwLT' 
'byVVTightIsolationMVA3oldDMwoLT' 'chargedIsoPtSum' 'decayModeFinding' 'decayModeFindingNewDMs' 'decayModeFindingOldDMs' 'neutralIsoPtSum' 
'puCorrPtSum' .
#endif

      // kinematic variables for PFJet associated to PFTau
      if (tb::isValidRef(it->pfJetRef())) {
        const reco::PFJetRef &jtr = it->pfJetRef();
        tau.jetPt  = jtr->pt();
        tau.jetEta = jtr->eta();
        tau.jetPhi = jtr->phi();
      }

      // NEW quantities
      tau.emFraction = it->emFraction();
      tau.maximumHCALPFClusterEt  = it->maximumHCALPFClusterEt();
      tau.ecalStripSumEOverPLead  = it->ecalStripSumEOverPLead();
      tau.bremsRecoveryEOverPLead = it->bremsRecoveryEOverPLead();
      tau.hcalTotOverPLead = it->hcalTotOverPLead();
      tau.hcalMaxOverPLead = it->hcalMaxOverPLead();
      tau.hcal3x3OverPLead = it->hcal3x3OverPLead();

      tau.etaetaMoment = it->etaetaMoment();
      tau.phiphiMoment = it->phiphiMoment();
      tau.phiphiMoment = it->etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      tau.vx = vertex.x();
      tau.vy = vertex.y();
      tau.vz = vertex.z();

      tau.zvertex = it->vz(); // distance from the primary vertex
      tau.mass    = it->p4().M();
      tau.ltsipt  = TMath::Abs(it->leadPFChargedHadrCandsignedSipt());

      // pat::Tau now has reference to the primary vertex also

      // add particle to the list
      list_->push_back(tau);
    }
    fnTau_ = list_->size();
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: "
                              << tauTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
