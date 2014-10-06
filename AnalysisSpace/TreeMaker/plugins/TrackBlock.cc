#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "Math/GenVector/VectorUtil.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/TrackBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

TrackBlock::TrackBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackSrc", edm::InputTag("generalTracks"))),
  bsTag_(iConfig.getUntrackedParameter<edm::InputTag>("offlineBeamSpot", edm::InputTag("offlineBeamSpot"))),
  trackToken_(consumes<reco::TrackCollection>(trackTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_))
{
}
void TrackBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Track>();
  tree->Branch("Track", "std::vector<vhtm::Track>", &list_);
  tree->Branch("nTrack", &fnTrack_, "fnTrack_/I");
}
void TrackBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  list_->clear();
  fnTrack_ = 0;

  // read the beam spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(bsToken_, beamSpot);

  // Read track collection
  edm::Handle<reco::TrackCollection> tracks;
  bool found = iEvent.getByToken(trackToken_, tracks);

  if (found && tracks.isValid()) {
    edm::LogInfo("TrackBlock") << "Total # of Tracks: " << tracks->size();
    reco::Track::TrackQuality quality = reco::Track::qualityByName("loose");
    for (auto it = tracks->begin(); it != tracks->end(); ++it) {
      const reco::Track &track = (*it);
      if (!track.quality(quality)) continue;

      if (fnTrack_ == kMaxTrack_) {
        edm::LogInfo("TrackBlock") << "Too many Tracks in the event, fnTrack = " 
                                   << fnTrack_;
        break;
      }
      vhtm::Track _tobj;
      _tobj.eta         = track.eta();
      _tobj.etaError    = track.etaError();
      _tobj.theta       = track.theta();
      _tobj.thetaError  = track.thetaError();
      _tobj.phi         = track.phi();
      _tobj.phiError    = track.phiError();
      _tobj.p           = track.p();
      _tobj.pt          = track.pt();
      _tobj.ptError     = track.ptError();
      _tobj.qoverp      = track.qoverp();
      _tobj.qoverpError = track.qoverpError();
      _tobj.charge      = track.charge();

      _tobj.nValidHits    = track.numberOfValidHits();
      _tobj.nLostHits     = track.numberOfLostHits();
      _tobj.validFraction = track.validFraction();

      const reco::HitPattern& hitp = track.hitPattern();
      _tobj.nValidTrackerHits            = hitp.numberOfValidTrackerHits();
      _tobj.nValidPixelHits              = hitp.numberOfValidPixelHits();
      _tobj.nValidStripHits              = hitp.numberOfValidStripHits();
      _tobj.trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();
      _tobj.pixelLayersWithMeasurement   = hitp.pixelLayersWithMeasurement();
      _tobj.stripLayersWithMeasurement   = hitp.stripLayersWithMeasurement();
 
      _tobj.dxy = track.dxy();
      _tobj.dz = track.dz();
      if (beamSpot.isValid()) {
        _tobj.dxy = track.dxy(beamSpot->position());
        _tobj.dz = track.dz(beamSpot->position());
      }
      _tobj.dxyError = track.dxyError();
      _tobj.dzError  = track.dzError();
     
      _tobj.chi2 = track.chi2();
      _tobj.ndof = track.ndof();
      _tobj.vx   = track.vx();
      _tobj.vy   = track.vy();
      _tobj.vz   = track.vz();

      list_->push_back(_tobj);
    }
    fnTrack_ = list_->size();
  }
  else {
    edm::LogError("TrackBlock") << "Error! Failed to get reco::Track collection, "
                                << trackTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackBlock);
