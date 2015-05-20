#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TTree.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/EventBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

EventBlock::EventBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  l1Tag_(iConfig.getUntrackedParameter<edm::InputTag>("l1Tag", edm::InputTag("gtDigis"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexTag", edm::InputTag("goodOfflinePrimaryVertices"))),
  pfTag_(iConfig.getUntrackedParameter<edm::InputTag>("pfTag", edm::InputTag("pfCands"))),
  selectedVertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("selectedVtxTag", edm::InputTag("selectedPrimaryVertices"))),
  puSummaryTag_(iConfig.getUntrackedParameter<edm::InputTag>("puSummaryTag", edm::InputTag("addPileupInfo"))),
  fixedGridRhoAllTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoAllTag", edm::InputTag("fixedGridRhoAll"))),
  fixedGridRhoFastjetAllTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllTag", edm::InputTag("fixedGridRhoFastjetAll"))),
  fixedGridRhoFastjetAllCaloTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetAllCaloTag", edm::InputTag("fixedGridRhoFastjetAllCalo"))),
  fixedGridRhoFastjetCentralCaloTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetCentralCaloTag", edm::InputTag("fixedGridRhoFastjetCentralCalo"))),
  fixedGridRhoFastjetCentralChargedPileUpTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetCentralChargedPileUpTag", edm::InputTag("fixedGridRhoFastjetCentralChargedPileUp"))),
  fixedGridRhoFastjetCentralNeutralTag_(iConfig.getUntrackedParameter<edm::InputTag>("fixedGridRhoFastjetCentralNeutralTag", edm::InputTag("fixedGridRhoFastjetCentralNeutral"))),
  vtxMinNDOF_(iConfig.getUntrackedParameter<unsigned int>("vertexMinimumNDOF", 4)),
  vtxMaxAbsZ_(iConfig.getUntrackedParameter<double>("vertexMaxAbsZ", 24.)),
  vtxMaxd0_(iConfig.getUntrackedParameter<double>("vertexMaxd0", 2.0)),
  numTracks_(iConfig.getUntrackedParameter<unsigned int>("numTracks", 10)),
  hpTrackThreshold_(iConfig.getUntrackedParameter<double>("hpTrackThreshold", 0.25)),
  l1Token_(consumes<L1GlobalTriggerReadoutRecord>(l1Tag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfTag_)),
  selectedVertexToken_(consumes<reco::VertexCollection>(selectedVertexTag_)),
  puSummaryToken_(consumes<std::vector<PileupSummaryInfo> >(puSummaryTag_)),
  fixedGridRhoAllToken_(consumes<double>(fixedGridRhoAllTag_)),
  fixedGridRhoFastjetAllToken_(consumes<double>(fixedGridRhoFastjetAllTag_)),
  fixedGridRhoFastjetAllCaloToken_(consumes<double>(fixedGridRhoFastjetAllCaloTag_)),
  fixedGridRhoFastjetCentralCaloToken_(consumes<double>(fixedGridRhoFastjetCentralCaloTag_)),
  fixedGridRhoFastjetCentralChargedPileUpToken_(consumes<double>(fixedGridRhoFastjetCentralChargedPileUpTag_)),
  fixedGridRhoFastjetCentralNeutralToken_(consumes<double>(fixedGridRhoFastjetCentralNeutralTag_)) 
{
}
EventBlock::~EventBlock() {
  delete nPU_;
  delete bunchCrossing_;
  delete trueNInt_;
  delete list_;
}
void EventBlock::beginJob() {

  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Event>();
  tree->Branch("Event", "std::vector<vhtm::Event>", &list_, 32000, 2);

  nPU_ = new std::vector<int>();
  tree->Branch("nPU", "std::vector<int>", &nPU_);

  bunchCrossing_ = new std::vector<int>();
  tree->Branch("bunchCrossing", "std::vector<int>", &bunchCrossing_);

  trueNInt_ = new std::vector<int>();
  tree->Branch("trueNInt", "std::vector<int>", &trueNInt_);
}
void EventBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector
  list_->clear();

  // Clear the independent vectors
  nPU_->clear();
  bunchCrossing_->clear();
  trueNInt_->clear();

  // Create Event Object
  vhtm::Event ev;
  ev.run   = iEvent.id().run();
  ev.event = iEvent.id().event();
  ev.lumis = iEvent.id().luminosityBlock();
  ev.bunch = iEvent.bunchCrossing();
  ev.orbit = iEvent.orbitNumber();

  double sec = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;
  ev.time   = sec+usec/conv;
  ev.isdata = iEvent.isRealData();
 
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  bool found = iEvent.getByToken(l1Token_, l1GtReadoutRecord);

  // Technical Trigger Part
  if (found && l1GtReadoutRecord.isValid()) {
    edm::LogInfo("EventBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: "
                               << l1Tag_;

    L1GtFdlWord fdlWord = l1GtReadoutRecord->gtFdlWord();
    if (fdlWord.physicsDeclared() == 1)
      ev.isPhysDeclared = true;

    // BPTX0
    if (l1GtReadoutRecord->technicalTriggerWord()[0])
      ev.isBPTX0 = true;

    // MinBias
    if (l1GtReadoutRecord->technicalTriggerWord()[40] || l1GtReadoutRecord->technicalTriggerWord()[41])
      ev.isBSCMinBias = true;

    // BeamHalo
    if ( (l1GtReadoutRecord->technicalTriggerWord()[36] || l1GtReadoutRecord->technicalTriggerWord()[37] ||
          l1GtReadoutRecord->technicalTriggerWord()[38] || l1GtReadoutRecord->technicalTriggerWord()[39]) ||
         ((l1GtReadoutRecord->technicalTriggerWord()[42] && !l1GtReadoutRecord->technicalTriggerWord()[43]) ||
          (l1GtReadoutRecord->technicalTriggerWord()[43] && !l1GtReadoutRecord->technicalTriggerWord()[42])) )
      ev.isBSCBeamHalo = true;
  }
  else {
    edm::LogError("EventBlock") << "Failed to get L1GlobalTriggerReadoutRecord for label: "
                                << l1Tag_;
  }

  // Good Primary Vertex Part
  edm::Handle<reco::VertexCollection> primaryVertices;
  found = iEvent.getByToken(vertexToken_, primaryVertices);

  if (found && primaryVertices.isValid()) {
    edm::LogInfo("EventBlock") << "Total # Primary Vertices: " << primaryVertices->size();
    for (const reco::Vertex& v: *primaryVertices) {
      if (!v.isFake() &&
           v.ndof() > vtxMinNDOF_ &&
	   std::fabs(v.z()) <= vtxMaxAbsZ_ &&
	   std::fabs(v.position().rho()) <= vtxMaxd0_)
      {
        ev.hasPrimaryVertex = true;
        break;
      }
    }
  }
  else {
    edm::LogError("EventBlock") << "Error! Failed to get VertexCollection for label: "
                                << vertexTag_;
  }
  edm::Handle<pat::PackedCandidateCollection> pfs;
  found = iEvent.getByToken(pfToken_, pfs);

  double sumPtPV = 0;
  int ntrk = 0, ntrkPV = 0;
  if (found && pfs.isValid()) {
    // now loop on pf candidates
    for (unsigned int i = 0; i < pfs->size(); ++i) {
      const pat::PackedCandidate& pf = (*pfs)[i];
      if (pf.charge() != 0) {
	++ntrk;
	if (pf.fromPV() == pat::PackedCandidate::PVAssoc::PVUsedInFit) {
	  ++ntrkPV;
	  sumPtPV += pf.pt();
	} 
      }
    }
  }
  ev.ntrk = ntrk;
  ev.ntrkPV = ntrkPV;
  ev.sumPtPV = sumPtPV;
#if 0
  // Scraping Events Part
  edm::Handle<reco::TrackCollection> tracks;
  found = iEvent.getByToken(trackToken_, tracks);

  if (found && tracks.isValid()) {
    edm::LogInfo("EventBlock") << "Total # Tracks: " << tracks->size();

    int numhighpurity = 0;
    double fraction = 1.;
    reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");
    if (tracks->size() > numTracks_) {
      for (const reco::Track& v: *tracks)
        if (v.quality(trackQuality)) numhighpurity++;
      fraction = static_cast<double>(numhighpurity)/static_cast<double>(tracks->size());
      if (fraction < hpTrackThreshold_)
        ev.isBeamScraping = true;
    }
  }
  else {
    edm::LogError("EventBlock") << "Error! Failed to get TrackCollection for label: "
                                << trackTag_;
  }
#endif
  // Access PU information
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
    found = iEvent.getByToken(puSummaryToken_, puInfo);
    if (found && puInfo.isValid()) {
      for (const PileupSummaryInfo& v: *puInfo) {
	ev.bunchCrossing.push_back(v.getBunchCrossing()); // in-time and out-of-time indices
	bunchCrossing_->push_back(v.getBunchCrossing());

	ev.nPU.push_back(v.getPU_NumInteractions());
	nPU_->push_back(v.getPU_NumInteractions());

	ev.trueNInt.push_back(v.getTrueNumInteractions());
	trueNInt_->push_back(v.getTrueNumInteractions());
      }
    }
    // More info about PU is here:
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
  }
  
  // Vertex Container
  edm::Handle<reco::VertexCollection> spVertices;
  found = iEvent.getByToken(selectedVertexToken_, spVertices);
  if (found) ev.nvtx = spVertices->size();

  // Rho
  edm::Handle<double> fixedGridRhoAll;
  iEvent.getByToken(fixedGridRhoAllToken_, fixedGridRhoAll);
  ev.fGridRhoAll = *fixedGridRhoAll;

  edm::Handle<double> fixedGridRhoFastjetAll;
  iEvent.getByToken(fixedGridRhoFastjetAllToken_, fixedGridRhoFastjetAll);
  ev.fGridRhoFastjetAll = *fixedGridRhoFastjetAll;

  edm::Handle<double> fixedGridRhoFastjetAllCalo;
  iEvent.getByToken(fixedGridRhoFastjetAllCaloToken_, fixedGridRhoFastjetAllCalo);
  ev.fGridRhoFastjetAllCalo = *fixedGridRhoFastjetAllCalo;

  edm::Handle<double> fixedGridRhoFastjetCentralCalo;
  iEvent.getByToken(fixedGridRhoFastjetCentralCaloToken_, fixedGridRhoFastjetCentralCalo);
  ev.fGridRhoFastjetCentralCalo = *fixedGridRhoFastjetCentralCalo;

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUp;
  iEvent.getByToken(fixedGridRhoFastjetCentralChargedPileUpToken_, fixedGridRhoFastjetCentralChargedPileUp);
  ev.fGridRhoFastjetCentralChargedPileUp = *fixedGridRhoFastjetCentralChargedPileUp;

  edm::Handle<double> fixedGridRhoFastjetCentralNeutral;
  iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken_, fixedGridRhoFastjetCentralNeutral);
  ev.fGridRhoFastjetCentralNeutral = *fixedGridRhoFastjetCentralNeutral;

  list_->push_back(ev);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventBlock);
