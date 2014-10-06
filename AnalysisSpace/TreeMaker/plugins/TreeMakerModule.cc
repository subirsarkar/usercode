#include <string>
#include <vector>
#include <cassert>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include "AnalysisSpace/TreeMaker/plugins/TreeMakerModule.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

TreeMakerModule::TreeMakerModule(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  createTree_(iConfig.getParameter<bool>("createTree"))
{
}
void TreeMakerModule::beginJob()
{
  if (!createTree_) return;
  edm::Service<TFileService> fs;
  fs->file().cd("/");
  TTree* tree = fs->make<TTree>("vhtree", "Physics Analysis Level TTree");
  assert(tree);
  fs->file().ls();
}
void TreeMakerModule::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get TTree pointer
  if (createTree_) return;
  TTree* tree = vhtm::Utility::getTree("vhtree");
  tree->Fill();
}
void TreeMakerModule::endJob() {
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TreeMakerModule);
