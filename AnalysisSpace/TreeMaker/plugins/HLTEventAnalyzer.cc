#include <cassert>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "AnalysisSpace/TreeMaker/plugins/HLTEventAnalyzer.h"

using namespace std;
using namespace edm;
 
HLTEventAnalyzer::HLTEventAnalyzer(const edm::ParameterSet& ps) : 
  processName_(ps.getParameter<std::string>("processName")),
  triggerName_(ps.getParameter<std::string>("triggerName")),
  triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
  triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent"))
{
  _srcTriggerResultsToken = consumes<edm::TriggerResults>(triggerResultsTag_);
  _srcTriggerEventToken = consumes<trigger::TriggerEvent>(triggerEventTag_);

  edm::LogInfo("HLTEventAnalyzer") << " configuration: \n"
				   << " ProcessName = " << processName_ << "\n"
				   << " TriggerName = " << triggerName_ << "\n"
				   << " TriggerResultsTag = " << triggerResultsTag_.encode() << "\n"
				   << " TriggerEventTag = " << triggerEventTag_.encode();
}
HLTEventAnalyzer::~HLTEventAnalyzer()
{
  edm::LogInfo("HLTEventAnalyzer") << " Destructor";
}
 
void HLTEventAnalyzer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      if (triggerName_ != "@") { // "@" means: analyze all triggers in config
	const unsigned int n(hltConfig_.size());
	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	if (triggerIndex >= n) {
	  edm::LogInfo("HLTEventAnalyzer") << "analyze:"
					   << " TriggerName " << triggerName_ 
					   << " not available in (new) config!";
	  edm::LogInfo("HLTEventAnalyzer") << "Available TriggerNames are: ";
	  hltConfig_.dump("Triggers");
	}
      }
      hltConfig_.dump("ProcessName");
      hltConfig_.dump("GlobalTag");
      hltConfig_.dump("TableName");
      hltConfig_.dump("Streams");
      hltConfig_.dump("Datasets");
      hltConfig_.dump("PrescaleTable");
      hltConfig_.dump("ProcessPSet");
    }
  } 
  else {
    edm::LogInfo("HLTEventAnalyzer") << "analyze:"
				     << " config extraction failure with process name "
				     << processName_;
  }
}
void HLTEventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::LogInfo("HLTEventAnalyzer")<< endl;
 
  // get event products
  iEvent.getByToken(_srcTriggerResultsToken, triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
     edm::LogInfo("HLTEventAnalyzer")<< "HLTEventAnalyzer::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  iEvent.getByToken(_srcTriggerEventToken, triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    edm::LogInfo("HLTEventAnalyzer")<< "HLTEventAnalyzer::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size() == hltConfig_.size());
   
  // analyze this event for the triggers requested
  if (triggerName_== "@") {
    const unsigned int n(hltConfig_.size());
    for (unsigned int i = 0; i != n; ++i) {
      analyzeTrigger(iEvent, iSetup, hltConfig_.triggerName(i));
    }
  } 
  else {
    analyzeTrigger(iEvent, iSetup, triggerName_);
  }
}
 
void HLTEventAnalyzer::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {
  using namespace reco;
  using namespace trigger;
 
  edm::LogInfo("HLTEventAnalyzer")<< endl;
 
  const unsigned int n(hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  assert(triggerIndex == iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
 
  // abort on invalid trigger name
  if (triggerIndex >= n) {
    edm::LogInfo("HLTEventAnalyzer") << "analyzeTrigger: path "
				     << triggerName << " - not found!" << endl;
    return;
  }
   
  const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerName));
  edm::LogInfo("HLTEventAnalyzer") << "analyzeTrigger: path "
				   << triggerName << " [" << triggerIndex << "] "
				   << "prescales L1T,HLT: " << prescales.first << "," << prescales.second;
 
  // modules on this trigger path
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  // Results from TriggerResults product
  edm::LogInfo("HLTEventAnalyzer") << " Trigger path status:"
				   << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
				   << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
				   << " Error =" << triggerResultsHandle_->error(triggerIndex);
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  edm::LogInfo("HLTEventAnalyzer") << " Last active module - label/type: "
				   << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
				   << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]";
  assert (moduleIndex < m);

  // Results from TriggerEvent product - Attention: must look only for
  // modules actually run in this path for this event!
  for (unsigned int j = 0; j <= moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      edm::LogInfo("HLTEventAnalyzer") << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType;
      const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
      const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      edm::LogInfo("HLTEventAnalyzer")<< "   " << n  << " accepted 'L3' objects found: " << endl;
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      for (size_type i = 0; i != n; ++i) {
	const TriggerObject& TO(TOC[KEYS[i]]);
	edm::LogInfo("HLTEventAnalyzer") << "  " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
					 << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass();
      }
    }
  }
}
