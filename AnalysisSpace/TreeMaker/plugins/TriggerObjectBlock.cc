#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "AnalysisSpace/TreeMaker/plugins/TriggerObjectBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

#include "TMath.h"
#include "TTree.h"
#include "TPRegexp.h"

// Constructor
TriggerObjectBlock::TriggerObjectBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  hltPathsOfInterest_(iConfig.getParameter<std::vector<std::string> >("hltPathsOfInterest")),
  hltPattern_(iConfig.getParameter<std::string>("hltPattern")),
  minTrigObjPt_(iConfig.getUntrackedParameter<double>("minTrigObjPt", 5.0)),
  hltTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults","","HLT"))),
  objectTag_(iConfig.getUntrackedParameter<edm::InputTag>("triggerObjectTag", edm::InputTag("selectedPatTrigger"))),
  hltToken_(consumes<edm::TriggerResults>(hltTag_)),
  objectToken_(consumes<pat::TriggerObjectStandAloneCollection>(objectTag_))
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
}
void TriggerObjectBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTriggerObject_ = 0;

  if (verbosity_) {
    std::cout << setiosflags(std::ios::fixed);
    std::cout << "Indx Eta Phi Pt Energy =Trigger path list=" << std::endl;
  }
 
  edm::Handle<edm::TriggerResults> triggerBits;
  bool found = iEvent.getByToken(hltToken_, triggerBits);
  if (found && triggerBits.isValid()) {
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    found = iEvent.getByToken(objectToken_, triggerObjects);
    if (found && triggerObjects.isValid()) {
      // Find the triggerNames and the matched paths
      const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
      for (const std::string& v: names.triggerNames()) 
	if (re_->Match(v)) matchedPathList_.push_back(v);

      for (pat::TriggerObjectStandAlone obj: *triggerObjects) {
	if (obj.pt() < minTrigObjPt_) continue;
	obj.unpackPathNames(names);

	std::map <std::string, unsigned int> pathInfoMap;
	for (const std::string& v: matchedPathList_) {
          int val = -1;
	  if (obj.hasPathName(v, true, true)) val = 3; 
	  else if (obj.hasPathName(v, false, true)) val = 2; 
	  else if (obj.hasPathName(v, true, false)) val = 1; 
	  else if (obj.hasPathName(v, false, false)) val = 0; 
	  if (val > -1) pathInfoMap.insert(std::pair<std::string, unsigned int>(v, val));
	}
	
	if (list_->size() == kMaxTriggerObject_) {
	  edm::LogInfo("TriggerObjectBlock") << "Too many Trigger Muons (HLT), fnTriggerObject = "
					     << list_->size();
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
    else {
      edm::LogError("TriggerObjectBlock") << "Failed to get TriggerObjects for label: "
					  << objectTag_;
    }
  }
  else {
    edm::LogError("TriggerObjectBlock") << "Failed to get TriggerResults for label: "
					<< hltTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
