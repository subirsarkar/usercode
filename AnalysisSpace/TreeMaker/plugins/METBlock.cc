#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/METBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

METBlock::METBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc", edm::InputTag("patMETs"))),
  corrMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("corrmetSrc", edm::InputTag("patMETsTypeIcorrected"))),
  mvaMETTag_(iConfig.getUntrackedParameter<edm::InputTag>("mvametSrc", edm::InputTag("patPFMetByMVA"))),
  pfMETToken_(consumes<pat::METCollection>(pfMETTag_)),
  corrMETToken_(consumes<pat::METCollection>(corrMETTag_)),
  mvaMETToken_(consumes<pat::METCollection>(mvaMETTag_))
{
}
void METBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  pfList_ = new std::vector<vhtm::MET>();
  tree->Branch("MET", "std::vector<vhtm::MET>", &pfList_, 32000, 2);
  tree->Branch("nMET", &fnPFMET_, "fnPFMET_/I");

#if 0
  corrList_ = new std::vector<vhtm::MET>();
  tree->Branch("corrMET", "std::vector<vhtm::MET>", &corrList_, 32000, 2);
  tree->Branch("corrnMET", &fnCorrMET_, "fnCorrMET_/I");

  mvaList_ = new std::vector<vhtm::MET>();
  tree->Branch("mvaMET", "std::vector<vhtm::MET>", &mvaList_, 32000, 2);
  tree->Branch("mvanMET", &fnMVAMET_, "fnMVAMET_/I");
#endif
}
void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  fillMET(iEvent, iSetup, pfList_, fnPFMET_, pfMETTag_, pfMETToken_);
  //fillMET(iEvent, iSetup, corrList_, fnCorrMET_, corrMETTag_, corrMETToken_);
  //fillMET(iEvent, iSetup, mvaList_, fnMVAMET_, mvaMETTag_, mvaMETToken_);
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<pat::METCollection>& token)
{
  // Reset the TClonesArray and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<pat::METCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const pat::MET& v: *metColl) {
      if (list->size() == kMaxMET_) {
        edm::LogInfo("METBlock") << "Too many PAT MET, nMET = " << list->size()
				 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      // fill in all the vectors
      mobj.met          = v.pt();
      mobj.metphi       = v.phi();
      mobj.sumet        = v.sumEt();
      mobj.metuncorr    = v.uncorrectedPt(pat::MET::uncorrALL);
      mobj.metphiuncorr = v.uncorrectedPhi(pat::MET::uncorrALL);
      mobj.sumetuncorr  = v.sumEt() - v.corSumEt(pat::MET::uncorrALL);

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get pat::MET collection for label: "
                              << iTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
