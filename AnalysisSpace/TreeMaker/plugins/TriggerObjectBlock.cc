#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

#include "AnalysisSpace/TreeMaker/plugins/TriggerObjectBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

#include "TMath.h"
#include "TTree.h"
#include "TPRegexp.h"

// Constructor
TriggerObjectBlock::TriggerObjectBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  hltTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults","","HLT"))),
  triggerEventTag_(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag", edm::InputTag("patTriggerEvent"))),
  hltPathsOfInterest_(iConfig.getParameter<std::vector<std::string> >("hltPathsOfInterest")),
  hltPattern_(iConfig.getParameter<std::string>("hltPattern")),
  minTrigObjPt_(iConfig.getUntrackedParameter<double>("minTrigObjPt", 8.0)),
  may10ReRecoData_(iConfig.getUntrackedParameter<bool>("May10ReRecoData", false)),
  triggerEventToken_(consumes<pat::TriggerEvent>(triggerEventTag_))
{
  std::cout << "hltPattern = " << std::endl
            << hltPattern_
            << std::endl;
  re_ = new TPMERegexp(hltPattern_, "xo");
}
TriggerObjectBlock::~TriggerObjectBlock() {
  if (re_) delete re_;
}
void TriggerObjectBlock::beginJob()
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::TriggerObject>();
  tree->Branch("TriggerObject", "std::vector<vhtm::TriggerObject>", &list_, 32000, 2);
  tree->Branch("nTriggerObject", &fnTriggerObject_, "fnTriggerObject_/I");

  firingFlag_ = (may10ReRecoData_) ? false : true;
}
void TriggerObjectBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig_.init(iRun, iSetup, hltTag_.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerObjectBlock") << "HLT config with process name "
                                       << hltTag_.process()
                                       << " successfully extracted";
    matchedPathList_.clear();
    auto list = hltConfig_.triggerNames();
    for (auto v: list) 
      if (re_->Match(v)) matchedPathList_.push_back(v);
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerObjectBlock") << "Error! HLT config extraction with process name "
                                        << hltTag_.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}
void TriggerObjectBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTriggerObject_ = 0;

  if (verbosity_) {
    std::cout << setiosflags(std::ios::fixed);
    std::cout << "Indx Eta Phi Pt Energy =Trigger path list=" << std::endl;
  }
 
  // trigger event
  edm::Handle<pat::TriggerEvent> triggerEvent;
  bool found = iEvent.getByToken(triggerEventToken_, triggerEvent);

  if (found && triggerEvent.isValid()) {
    // get the trigger objects corresponding to the used matching (HLT muons) and
    // loop over selected trigger objects
    pat::TriggerObjectRefVector objectList(triggerEvent->objectRefs());
    for (auto it = objectList.begin(); it != objectList.end(); ++it) {
      const pat::TriggerObject& obj = (**it);
      if (obj.pt() < minTrigObjPt_) continue;
      
      std::map <std::string, unsigned int> pathInfoMap;
      for (auto jt: matchedPathList_) {
	std::string name = jt;
	if (!triggerEvent->objectInPath(*it, name)) continue;
	
	bool matched = true;
	// Get the filters and access the L3 filter (needed for May10ReReco data)
	if (may10ReRecoData_) {
	  matched = false;
	  pat::TriggerFilterRefVector filters(triggerEvent->pathFilters(name, firingFlag_));
	  if (filters.empty()) continue;
	  pat::TriggerFilterRef lastFilter(filters.at(filters.size() - 1));
	  if (triggerEvent->objectInFilter(*it, lastFilter->label())) matched = true;
	}
	if (matched) {
	  unsigned int val = (triggerEvent->path(name)->wasRun() && triggerEvent->path(name)->wasAccept()) ? 1 : 0;
	  pathInfoMap.insert(std::pair<std::string, unsigned int>(name, val));
	}
      }
      if (pathInfoMap.size() > 0) {
	if (fnTriggerObject_ == kMaxTriggerObject_) {
	  edm::LogInfo("TriggerObjectBlock") << "Too many Trigger Muons (HLT), fnTriggerObject = "
					     << fnTriggerObject_;
	  break;
	}
	vhtm::TriggerObject _tobj;
	_tobj.eta      = obj.eta();
	_tobj.phi      = obj.phi();
	_tobj.pt       = obj.pt();
	_tobj.energy   = obj.energy();
	_tobj.pathList = pathInfoMap;
	
	if (verbosity_) {
	  std::cout << std::setprecision(2);
	  std::cout << std::setw(4) << fnTriggerObject_++
		    << std::setw(8) << _tobj.eta
		    << std::setw(8) << _tobj.phi
		    << std::setw(8) << _tobj.pt
		    << std::setw(8) << _tobj.energy
		    << std::endl;
	  for (auto jt: _tobj.pathList)
	    std::cout << "\t\t\t\t\t" << jt.first << " " << jt.second << std::endl;
	}
	list_->push_back(_tobj);
      }
      fnTriggerObject_ = list_->size();
    }
    if (verbosity_) std::cout << " # of Trigger Objects " << fnTriggerObject_ << std::endl;
  }
  else {
    edm::LogError("TriggerObjectBlock") << "Failed to get TriggerResults for label: "
					<< triggerEventTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
