#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/VertexBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

VertexBlock::VertexBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc", edm::InputTag("goodOfflinePrimaryVertices"))),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_))
{
}
void VertexBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Vertex>();
  tree->Branch("Vertex", "std::vector<vhtm::Vertex>", &list_, 32000, 2);
  tree->Branch("nVertex", &fnVertex_, "fnVertex_/I");
}
void VertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnVertex_ = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  bool found = iEvent.getByToken(vertexToken_, primaryVertices);

  if (found && primaryVertices.isValid()) {
    edm::LogInfo("VertexBlock") << "Total # of Primary Vertices: " 
                                << primaryVertices->size();

    for (auto it = primaryVertices->begin(); it != primaryVertices->end(); ++it) {
      if (fnVertex_ == kMaxVertex_) {
        edm::LogInfo("VertexBlock") << "Too many Vertex, fnVertex = " 
                                    << fnVertex_;
        break;
      }
      vhtm::Vertex vertex;
      vertex.x          = it->x();
      vertex.y          = it->y();
      vertex.z          = it->z();
      vertex.xErr       = it->xError();
      vertex.yErr       = it->yError();
      vertex.zErr       = it->zError();
      vertex.rho        = it->position().rho();
      vertex.chi2       = it->chi2();
      vertex.ndf        = it->ndof();
      vertex.ntracks    = static_cast<int>(it->tracksSize());
      vertex.ntracksw05 = it->nTracks(0.5); // number of tracks in the vertex with weight above 0.5
      vertex.isfake     = it->isFake();
      vertex.isvalid    = it->isValid();
      vertex.sumPt      = it->p4().pt();

      list_->push_back(vertex);
    }
    fnVertex_ = list_->size();
  }
  else {
    edm::LogError("VertexBlock") << "Error! Failed to get VertexCollection for label: "
                                 << vertexTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexBlock);
