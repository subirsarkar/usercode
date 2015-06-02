#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/METBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

METBlock::METBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc", edm::InputTag("patMETs"))),
  mvaMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvametSrc", edm::InputTag("pfMVAMEt"))),
  pfMETToken_(consumes<pat::METCollection>(pfMETTag_)),
  mvaMETToken_(consumes<reco::PFMETCollection>(mvaMETTag_))
{
}
METBlock::~METBlock() {
  delete pfList_;
}
void METBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  pfList_ = new std::vector<vhtm::MET>();
  tree->Branch("MET", "std::vector<vhtm::MET>", &pfList_, 32000, 2);
  tree->Branch("nMET", &fnPFMET_, "fnPFMET_/I");

  mvaList_ = new std::vector<vhtm::MET>();
  tree->Branch("mvaMET", "std::vector<vhtm::MET>", &mvaList_, 32000, 2);
  tree->Branch("mvanMET", &fnMVAMET_, "fnMVAMET_/I");
}
void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  fillMET(iEvent, iSetup, pfList_, fnPFMET_, pfMETTag_, pfMETToken_);
  fillMET(iEvent, iSetup, mvaList_, fnMVAMET_, mvaMETTag_, mvaMETToken_);
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<pat::METCollection>& token)
{
  // Reset the vector and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<pat::METCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const pat::MET& v: *metColl) {
      if (list->size() == kMaxMET_) {
        edm::LogInfo("METBlock") << "Too many pat::MET, nMET = " << list->size()
				 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      mobj.met          = v.pt();
      mobj.metphi       = v.phi();
      mobj.sumet        = v.sumEt();
      mobj.metuncorr    = v.uncorrectedPt(pat::MET::uncorrALL);
      mobj.metphiuncorr = v.uncorrectedPhi(pat::MET::uncorrALL);
      mobj.sumetuncorr  = v.sumEt() - v.corSumEt(pat::MET::uncorrALL);
      mobj.metJESUp     = v.shiftedPt(pat::MET::JetEnUp);
      mobj.metJESDn     = v.shiftedPt(pat::MET::JetEnDown);

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get pat::MET collection for label: "
                              << iTag;
  }
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<reco::PFMETCollection>& token)
{
  // Reset the vector and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<reco::PFMETCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const reco::PFMET& v: *metColl) {
      if (list->size() == kMaxMET_) {
	edm::LogInfo("METBlock") << "Too many PFMET, nMET = " << list->size()
                                 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      mobj.met    = v.pt();
      mobj.metphi = v.phi();
      mobj.sumet  = v.sumEt();

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get reco::PFMET collection for label: "
                              << iTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
