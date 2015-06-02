#ifndef __AnalysisSpace_TreeMaker_PhysicsObjects_h
#define __AnalysisSpace_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>

#include "TObject.h"

namespace vhtm {
  class Candidate: public TObject {
  public:
    Candidate();
    Candidate(float pt, float eta, float phi);
    virtual ~Candidate() {}

    float pt;
    float eta;
    float phi;

    ClassDef(Candidate, 1)
  };
  class PackedPFCandidate: public TObject {
  public:
    PackedPFCandidate();
    virtual ~PackedPFCandidate() {}

    float pt;
    float eta;
    float phi;
    float energy;
    
    int pdgId;
    int charge;
    
    double vx;
    double vy;
    double vz;
   
    int fromPV;
    // w.r.t PV
    float dxy;
    float dz;
    float dxyError;
    float dzError;

    std::map<std::string, std::vector<double> > isolationMap;

    ClassDef(PackedPFCandidate, 1)
  };
  class Event: public TObject {
  public:
    Event();
    virtual ~Event() {}
  
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
    bool hasPrimaryVertex;
    int ntrk;
    int ntrkPV;
    int sumPtPV;
    //bool isBeamScraping;
    bool passHBHENoiseFilter;
  
    double fGridRhoAll;
    double fGridRhoFastjetAll;
    double fGridRhoFastjetAllCalo;
    double fGridRhoFastjetCentralCalo;
    double fGridRhoFastjetCentralChargedPileUp;
    double fGridRhoFastjetCentralNeutral;
    int nvtx; 

    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;
  
    ClassDef(Event, 1)
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    virtual ~GenEvent() {}
  
    unsigned int processID;
    double ptHat;
    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent, 1)
  };
  class Electron: public TObject {
  public:
    Electron();
    virtual ~Electron() {}
    float eta;
    float phi;
    float pt;
    bool ecalDriven;
    bool hasGsfTrack;
    float energy;
    float caloEnergy;
    float charge;
    float eOverPOut;

    float trackP;
    float trackPt;
    double trkD0;
    double trkDz;
    float normalizedChi2;
    int pixHits;
    int trkHits;
    int nValidHits;
    int missingHits;
    int nLayers;

    double dxyPV;
    double dzPV;

    double vtxDist3D;
    int vtxIndex;
    double vtxDistZ;

    // ID variables
    float hoe;
    float eop;
    float sigmaIEtaIEta;
    float sigmaIPhiIPhi;
    float deltaPhiTrkSC;
    float deltaEtaTrkSC;
    float deltaEtaCalo;
    float r9;
    float e1x5;
    float e5x5;
  
    // SC associated with electron
    float scEn;
    float scEta;
    float scPhi;
    float scRawEnergy;
    float scEtaWidth;
    float scPhiWidth;
    float scPreshowerEnergy;

    // Vertex association variables
    float pfRelIso;

    // PFlow isolation variable
    float chargedHadronIso;
    float neutralHadronIso;
    float photonIso;

    // PFlow isolation variable
    float sumChargedParticlePt;
    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;
 
    // IP against PV
    double dB;
    double edB;
    double dB3D;
    double edB3D;

    // Bremsstrahlung
    int nBrems;
    float fbrem;
  
    bool hasMatchedConv;
    bool mvaPreselection;
    bool isTriggerElectron;
    int fidFlag;

    // Vertex
    float vx;
    float vy;
    float vz;

    std::map<std::string, float> idmap;
    float mvaId; 
    int selbit;

    std::map<std::string, std::vector<double> > isolationMap;

    ClassDef(Electron, 1)
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    virtual ~GenParticle() {}
  
    float eta;
    float phi;
    float p;
    float px;
    float py;
    float pz;
    float pt;
    float energy;
    int pdgId;
    float vx;
    float vy;
    float vz;
    int status;
    float charge;
    int numDaught;
    int numMother;
    int motherIndex;
    std::vector<int> motherIndices;
    std::vector<int> daughtIndices;
  
    ClassDef(GenParticle, 1)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    virtual ~GenJet() {}
  
    float eta;
    float phi;
    float p;
    float pt;
    float energy;
    float emf;
    float hadf;
  
    ClassDef(GenJet, 1)
  };
  class MET: public TObject {
  public:
    MET();
    virtual ~MET() {}
  
    float met;
    float metphi;
    float sumet;
    float metuncorr;
    float metphiuncorr;
    float sumetuncorr;
    float metJESUp;
    float metJESDn;
  
    ClassDef(MET, 1)
  };
  class Tau: public TObject {
  public:
    Tau();
    virtual ~Tau() {}
  
    float eta;
    float phi;
    float pt;
    float energy;
    float charge;
    float mass;
  
    double dxyPV;
    double dzPV;
    int vtxIndex;
    double vtxDxy;
    double vtxDz;

    // Leading particle pT
    float leadChargedParticlePt;
    float leadNeutralParticlePt;
    float leadParticlePt;

    std::vector<vhtm::Candidate> sigChHadList;
    std::vector<vhtm::Candidate> sigNeHadList;
    std::vector<vhtm::Candidate> sigGammaList;
    std::vector<vhtm::Candidate> isoChHadList;
    std::vector<vhtm::Candidate> isoNeHadList;
    std::vector<vhtm::Candidate> isoGammaList;

    float ptSumChargedHadronsIsoCone;
    float ptSumNeutralHadronsIsoCone;
    float ptSumPhotonsIsoCone;

    // tau id. discriminators
    float decayModeFinding;
    float decayModeFindingNewDMs;
    float decayModeFindingOldDMs;
    
    // discriminators against electrons/muons
    float againstMuonLoose;
    float againstMuonMedium;
    float againstMuonTight;
    
    float againstMuonLoose3;
    float againstMuonTight3;

    float againstElectronLoose;
    float againstElectronMedium;
    float againstElectronTight;
     //float againstElectronMVA;
     
    float againstElectronLooseMVA5;
    float againstElectronMediumMVA5;
    float againstElectronTightMVA5;
    
    float byLooseCombinedIsolationDeltaBetaCorr3Hits;
    float byMediumCombinedIsolationDeltaBetaCorr3Hits;
    float byTightCombinedIsolationDeltaBetaCorr3Hits;
    float byCombinedIsolationDeltaBetaCorrRaw3Hits;
    float chargedIsoPtSum;
    float neutralIsoPtSum;
    float puCorrPtSum;
    
    // kinematic variables for PFJet associated to PFTau
    float jetPt;
    float jetEta;
    float jetPhi;
    float emFraction;
    float vx;
    float vy;
    float vz;
    
    float zvertex;
    float dxySig;
    int selbit;
    
    ClassDef(Tau, 1)
  };
  class Muon: public TObject {
  public:
    Muon();
    virtual ~Muon() {}
    bool isGlobalMuon;
    bool isTrackerMuon;
    bool isPFMuon;
    bool isGhostCleaned;

    float eta;
    float phi;
    float pt;
    float p;
    float energy;
    float charge;
    bool passID;

    int muonBestTrackType;  
    double trkD0;
    double trkDz;
    float normChi2;

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
    int numMuonStations;

    float trkIso;
    float ecalIso;
    float hcalIso;
    float hoIso;

    float sumChargedParticlePtR03;
    float sumChargedHadronPtR03;
    float sumNeutralHadronEtR03;
    float sumPhotonEtR03;
    float sumPUPtR03;
    float pfRelIso03;

    float sumChargedParticlePt;
    float sumChargedHadronPt;
    float sumNeutralHadronEt;
    float sumPhotonEt;
    float sumPUPt;
    float pfRelIso04;

    double dB; // PV2D
    double edB; 
    double dB3D;
    double edB3D;
  
    // UW Recommendation
    bool isGlobalMuonPromptTight;
    bool isAllArbitrated;
    int nChambers;
    int nMatches;
    int nMatchedStations;
    unsigned int stationMask;
    unsigned int stationGapMaskDistance;
    unsigned int stationGapMaskPull;
    bool muonID;

    float vx;
    float vy;
    float vz;

    std::map<std::string, std::vector<double> > isolationMap;
    int nSegments;

    int selbit;

    ClassDef(Muon, 1)
  };
  class Jet: public TObject {
  public:
    Jet();
    virtual ~Jet() {}
    float eta;
    float phi;
    float pt;
    float pt_raw;
    float energy;
    float energy_raw;
    float jecUnc;
    float resJEC;
    int partonFlavour;

    float chargedEmEnergyFraction;
    float chargedHadronEnergyFraction;
    float chargedMuEnergyFraction;
    float electronEnergyFraction;
    float muonEnergyFraction;
    float neutralEmEnergyFraction;
    float neutralHadronEnergyFraction;
    float photonEnergyFraction;
    int chargedHadronMultiplicity;
    int chargedMultiplicity;
    int electronMultiplicity;
    int muonMultiplicity;
    int neutralHadronMultiplicity;
    int neutralMultiplicity;
    int photonMultiplicity;
    int nConstituents;

    float combinedSecondaryVertexBTag;
    float combinedInclusiveSecondaryVertexV2BJetTags;
    float combinedInclusiveSecondaryVertexBTag;

    std::map<std::string, float> discrimap;

    int passLooseID;
    int passTightID;
  
    float jpumva;
    int selbit;

    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    virtual ~Vertex() {}
  
    float x;
    float y;
    float z;
    float xErr;
    float yErr;
    float zErr;
    float rho;
    float chi2;
    float ndf;
    int ntracks;
    bool isfake;
    bool isvalid;

    int selbit;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    virtual ~GenMET() {}
  
    float met;
    float metphi;
    float sumet;
  
    ClassDef(GenMET, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    virtual ~TriggerObject() {}
  
    float pt;
    float eta;
    float phi;
    float energy;
  
    std::map<std::string, unsigned int> pathList;
  
    ClassDef(TriggerObject, 1)
  };
  class Photon : public TObject {
  public:
    Photon();
    virtual ~Photon() {}
    
    float et;
    float eta;
    float clusterEta;
    float phi;
    float clusterPhi;
    float energy;
    float theta; 
    float vx;
    float vy;
    float vz;
    
    float scEnergy;
    float scEta;
    float scPhi;
    float scSize;
    float scEtaWidth;
    float scPhiWidth;
    float scEt;
    float scRawEnergy;
    float scx;
    float scy;
    float scz; 

    float isoEcalRecHit03;
    float isoHcalRecHit03;
    float isoSolidTrkCone03;
    float isoHollowTrkCone03;
    int nTrkSolidCone03;
    int nTrkHollowCone03;
    
    float isoEcalRecHit04;
    float isoHcalRecHit04;
    float isoSolidTrkCone04;
    float isoHollowTrkCone04;
    int nTrkSolidCone04;
    int nTrkHollowCone04;
    
    bool isEB;
    bool isEE; 
    bool isEBGap;
    bool isEEGap;
    bool isEBEEGap;
    int fidFlag;
    
    bool hasPixelSeed;
    bool passElectronVeto;

    float chargedHadIso;
    float neutralHadIso;
    float photonIso;
    float puChargedHadIso;
    
    float r9;
    float hoe;
    float sigmaEtaEta;
    float sigmaIEtaIEta;
    float e1x5;
    float e2x5; 
    float e3x3;
    float e5x5;
    float r1x5;
    float r2x5;
    float maxEnergyXtal;
    std::map<std::string, float> idmap;
    
    bool hasConversionTracks;
    int nTracks;
    bool isConverted;
    float pairInvMass;
    float pairCotThetaSeparation;
    float pairPx;
    float pairPy;
    float pairPz;
    float conv_vx;
    float conv_vy;
    float conv_vz;
    float eovp;
    float zpv;
    float distOfMinApproach;
    float dPhiTracksAtVtx;
    float dPhiTracksAtEcal;
    float dEtaTracksAtEcal;  
    
    int selbit;
  
    ClassDef(Photon, 1)
  };
}
#endif
