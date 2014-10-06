#ifndef __AnalysisSpace_TreeMaker_PhysicsObjects_h
#define __AnalysisSpace_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>

#include "TObject.h"

namespace vhtm {
  class Event: public TObject {
  public:
    Event();
    ~Event() {}
  
    unsigned int run;
    unsigned int event;
    unsigned int lumis;
    unsigned int bunch;
    unsigned int orbit;
          double time;
            bool isdata;
  
    bool isPhysDeclared;
    bool isBPTX0;
    bool isBSCMinBias;
    bool isBSCBeamHalo;
    bool isPrimaryVertex;
    bool isBeamScraping;
    bool passHBHENoiseFilter;
  
    double rho;
    double rhoNeutral;
    int nvtx; 

    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;
  
    ClassDef(Event,1)
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    ~GenEvent() {}
  
    unsigned int processID;
    double ptHat;
    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent,1)
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}
    double eta;
    double phi;
    double pt;
      bool ecalDriven;
      bool hasGsfTrack;
    float trackPt;
    double trackPtError;
    double energy;
    double caloEnergy;
    double caloEnergyError;
       int charge;
       int pixHits;
       int trkHits;
       int nValidHits;
    double trkD0;
    double trkD0Error;
    double trkDz;
    double trkDzError;
  
    // ID variables
    float hoe;
    float hoeDepth1;
    float eop;
    float sigmaEtaEta;
    float sigmaIEtaIEta;
    float deltaPhiTrkSC;
    float deltaEtaTrkSC;
      int classif;
    float e1x5overe5x5;
    float e2x5overe5x5;
  
    // Iso variables
    float isoEcal03;
    float isoHcal03;
    float isoTrk03;
    float isoEcal04;
    float isoHcal04;
    float isoTrk04;
    float isoRel03;
    float isoRel04;
  
    // Vertex
    double vx;
    double vy;
    double vz;

    // SC associated with electron
    double scEn;
    double scEta;
    double scPhi;
    double scET;
    double scRawEnergy;
  
    // Vertex association variables
    double dxyPV;
    double dzPV;
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
    double relIso;
#if 1
    double pfRelIso;
#endif  
    // PFlow isolation variable
    float chargedHadronIso;
    float neutralHadronIso;
    float photonIso;

    float sumChargedHadronPt;
    float sumPUPt;

    int missingHits;

    double dB;
    double edB;
  
    double dB3d;
    double edB3d;
  
    double scE1E9;
    double scS4S1;
    double sckOutOfTime;
    double scEcalIso;
    double scHEEPEcalIso;
    double scHEEPTrkIso;
  
       int nBrems;
     float fbrem;
  
     float dist_vec;
     float dCotTheta;
      bool hasMatchedConv;
#if 1    
     float mva;
     float mvaPOGTrig;
     float mvaPOGNonTrig;
#endif
     bool mvaPreselection;
     bool isTriggerElectron;
#if 0
     float isoMVA;
     float pfRelIso03v1;
     float pfRelIso03v2;
     float pfRelIsoDB03v1;
     float pfRelIsoDB03v2;
     float pfRelIsoDB03v3;

     float pfRelIso04v1;
     float pfRelIso04v2;
     float pfRelIsoDB04v1;
     float pfRelIsoDB04v2;
     float pfRelIsoDB04v3;

     float pfRelIso03;
     float pfRelIso04;
     float pfRelIsoDB03;
     float pfRelIsoDB04;
#endif   

       int selbit;
       int fidFlag;

    ClassDef(Electron, 3)
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
  
    double eta;
    double phi;
    double p;
    double px;
    double py;
    double pz;
    double pt;
    double energy;
       int pdgId;
    double vx;
    double vy;
    double vz;
       int status;
    double charge;
       int numDaught;
       int numMother;
       int motherIndex;
    std::vector<int> motherIndices;
    std::vector<int> daughtIndices;
  
    ClassDef(GenParticle,1)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    float eta;
    float phi;
    float p;
    float pt;
    float energy;
    float emf;
    float hadf;
  
    ClassDef(GenJet,1)
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    double met;
    double metphi;
    double sumet;
    double metuncorr;
    double metphiuncorr;
    double sumetuncorr;
  
    ClassDef(MET,1)
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    enum {
      kMaxPFChargedCand = 40,
      kMaxPFNeutralCand = 20
    };
    double eta;
    double phi;
    double pt;
    double energy;
       int charge;
    double mass;
  
    double leadTrkPt;
    double leadTrkPtError;
    double leadTrkEta;
    double leadTrkPhi;
    double leadTrkCharge;
    double leadTrkD0;
    double leadTrkD0Error;
    double leadTrkDz;
    double leadTrkDzError;

    double dxyPV;
    double dzPV;
       int vtxIndex;
    double vtxDxy;
    double vtxDz;

    // Leading particle pT
    double leadChargedParticlePt;
    double leadNeutralParticlePt;
    double leadParticlePt;
  
    // Number of charged/neutral candidates and photons in different cones
       int numChargedHadronsSignalCone;
       int numNeutralHadronsSignalCone;
       int numPhotonsSignalCone;
       int numParticlesSignalCone;
       int numChargedHadronsIsoCone;
       int numNeutralHadronsIsoCone;
       int numPhotonsIsoCone;
       int numParticlesIsoCone;
    double ptSumChargedHadronsIsoCone;
    double ptSumNeutralHadronsIsoCone;
    double ptSumPhotonsIsoCone;
    
    double sigChHadCandPt[kMaxPFChargedCand];
    double sigChHadCandEta[kMaxPFChargedCand];
    double sigChHadCandPhi[kMaxPFChargedCand];
    double sigNeHadCandPt[kMaxPFNeutralCand];
    double sigNeHadCandEta[kMaxPFNeutralCand];
    double sigNeHadCandPhi[kMaxPFNeutralCand];
    double sigGammaCandPt[kMaxPFNeutralCand];
    double sigGammaCandEta[kMaxPFNeutralCand];
    double sigGammaCandPhi[kMaxPFNeutralCand];

    double isoChHadCandPt[kMaxPFChargedCand];
    double isoChHadCandEta[kMaxPFChargedCand];
    double isoChHadCandPhi[kMaxPFChargedCand];
    double isoNeHadCandPt[kMaxPFNeutralCand];
    double isoNeHadCandEta[kMaxPFNeutralCand];
    double isoNeHadCandPhi[kMaxPFNeutralCand];
    double isoGammaCandPt[kMaxPFNeutralCand];
    double isoGammaCandEta[kMaxPFNeutralCand];
    double isoGammaCandPhi[kMaxPFNeutralCand];

     // tau id. discriminators
     float decayModeFinding;
     float decayModeFindingNewDMs;
     float decayModeFindingOldDMs;

     // discriminators against electrons/muons
     float againstMuonLoose;
     float againstMuonMedium;
     float againstMuonTight;

     float againstMuonLoose2;
     float againstMuonMedium2;
     float againstMuonTight2;

     float againstMuonLoose3;
     float againstMuonTight3;

     float againstElectronLoose;
     float againstElectronMedium;
     float againstElectronTight;
     float againstElectronMVA;
  
     float againstElectronLooseMVA5;
     float againstElectronMediumMVA5;
     float againstElectronTightMVA5;

     float againstElectronVLooseMVA2;
     float againstElectronLooseMVA2;
     float againstElectronMediumMVA2;
     float againstElectronTightMVA2;

     float againstElectronLooseMVA3;
     float againstElectronMediumMVA3;
     float againstElectronTightMVA3;
     float againstElectronVTightMVA3;

     // Obsolete
     float pfElectronMVA;

     float byVLooseIsolation;
     float byLooseIsolation;
     float byMediumIsolation;
     float byTightIsolation;
  
     float byLooseIsolationMVA;
     float byMediumIsolationMVA;
     float byTightIsolationMVA;

     float byLooseIsolationMVA2;
     float byMediumIsolationMVA2;
     float byTightIsolationMVA2;

     float byVLooseCombinedIsolationDeltaBetaCorr;
     float byLooseCombinedIsolationDeltaBetaCorr;
     float byMediumCombinedIsolationDeltaBetaCorr;
     float byTightCombinedIsolationDeltaBetaCorr;
     float byVLooseIsolationDeltaBetaCorr;
     float byLooseIsolationDeltaBetaCorr;
     float byMediumIsolationDeltaBetaCorr;
     float byTightIsolationDeltaBetaCorr;
  
     float byLooseCombinedIsolationDeltaBetaCorr3Hits;
     float byMediumCombinedIsolationDeltaBetaCorr3Hits;
     float byTightCombinedIsolationDeltaBetaCorr3Hits;
     float byCombinedIsolationDeltaBetaCorrRaw3Hits;
     float chargedIsoPtSum;
     float neutralIsoPtSum;
     float puCorrPtSum;

     // kinematic variables for PFJet associated to PFTau
    double jetPt;
    double jetEta;
    double jetPhi;
  
     float emFraction;
     float maximumHCALPFClusterEt;
     float ecalStripSumEOverPLead;
     float bremsRecoveryEOverPLead;
     float hcalTotOverPLead;
     float hcalMaxOverPLead;
     float hcal3x3OverPLead;
  
     float etaetaMoment;
     float phiphiMoment;
     float etaphiMoment;
  
    double vx;
    double vy;
    double vz;
  
    double zvertex;
    double ltsipt;
  
       int selbit;
  
    ClassDef(Tau, 3)
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
      bool isTrackerMuon;
      bool isPFMuon;
    double eta;
    double phi;
    double pt;
    double ptError;
    double p;
    double energy;
       int charge;
    double trkD0;
    double trkD0Error;
    double trkDz;
    double trkDzError;
    double globalChi2;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double hoIso;
    double relIso;
    double pfChargedIsoR03;
    double pfChargedIsoR04;
       int passID;
    double dxyPV;
    double dzPV;
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
       int pixHits;
       int trkHits;
       int muoHits;
       int matches;
       int trackerLayersWithMeasurement;
#if 0
    double pfRelIso;
#endif
    double vx;
    double vy;
    double vz;

    double dB; // PV2D
    double edB;

    double dB3d; // PV3D
    double edB3d;
  
    // UW Recommendation
      bool isGlobalMuonPromptTight;
      bool isAllArbitrated;
       int nChambers;
       int nMatches;
       int nMatchedStations;
       unsigned int stationMask;
       unsigned int stationGapMaskDistance;
       unsigned int stationGapMaskPull;
       int muonID;
#if 0     
       float idMVA;
       float isoRingsMVA;
       float isoRingsRadMVA;
       float idIsoCombMVA;

       float pfRelIso03v1;
       float pfRelIso03v2;
       float pfRelIsoDB03v1;
       float pfRelIsoDB03v2;
       float pfRelIso04v1;
       float pfRelIso04v2;
       float pfRelIsoDB04v1;
       float pfRelIsoDB04v2;
#endif
       int selbit;

     ClassDef(Muon, 3)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
       int partonFlavour;
#if 0
     float puIdMVA;
       int puIdFlag;
       int puIdBits;
#endif
    double chargedEmEnergyFraction;
    double chargedHadronEnergyFraction;
    double chargedMuEnergyFraction;
    double electronEnergyFraction;
    double muonEnergyFraction;
    double neutralEmEnergyFraction;
    double neutralHadronEnergyFraction;
    double photonEnergyFraction;
       int chargedHadronMultiplicity;
       int chargedMultiplicity;
       int electronMultiplicity;
       int muonMultiplicity;
       int neutralHadronMultiplicity;
       int neutralMultiplicity;
       int photonMultiplicity;
       int nConstituents;
    double trackCountingHighEffBTag;
    double trackCountingHighPurBTag;
    double simpleSecondaryVertexHighEffBTag;
    double simpleSecondaryVertexHighPurBTag;
    double jetProbabilityBTag;
    double jetBProbabilityBTag;
    double combinedSecondaryVertexBTag;
    double combinedSecondaryVertexMVABTag;
    double combinedInclusiveSecondaryVertexBTag;
    double combinedMVABTag;
       int passLooseID;
       int passTightID;
  
       int selbit;

    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    ~Vertex() {}
  
    double x;
    double y;
    double z;
    double xErr;
    double yErr;
    double zErr;
    double rho;
    double chi2;
    double ndf;
       int ntracks;
       int ntracksw05;
      bool isfake;
      bool isvalid;
    double sumPt; // vector sum
  
       int selbit;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class Track: public TObject {
  public:
    Track();
    ~Track() {}
  
    double eta;
    double etaError;
    double theta;
    double thetaError;
    double phi;
    double phiError;
    double p;
    double pt;
    double ptError;
    double qoverp;
    double qoverpError;
     float charge;

       int nValidHits;
       int nLostHits;
    double validFraction;
       int nValidTrackerHits;
       int nValidPixelHits;
       int nValidStripHits;
       int trackerLayersWithMeasurement;
       int pixelLayersWithMeasurement;
       int stripLayersWithMeasurement;
  
    double dxy;
    double dxyError;
    double dz;
    double dzError;
  
    double chi2;
       int ndof;
    double vx;
    double vy;
    double vz;
  
       int selbit;
         
    ClassDef(Track, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    ~TriggerObject() {}
  
    double energy;
    double pt;
    double eta;
    double phi;
  
    std::map< std::string, unsigned int > pathList;
  
    ClassDef(TriggerObject, 1)
  };
}
#endif
