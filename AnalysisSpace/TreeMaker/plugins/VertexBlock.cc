#include <iostream>
#include <iomanip>
#include <sstream>

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
{}
VertexBlock::~VertexBlock() {
  delete list_;
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

    if (verbosity_) {
      edm::LogInfo("VertexBlock") << std::setprecision(2); 
      edm::LogInfo("VertexBlock") << "   indx      x      y      z    rho     chi2      ndf ntracks" << std::endl;
    }
    for (const reco::Vertex& v: *primaryVertices) {
      if (list_->size() == kMaxVertex_) {
        edm::LogInfo("VertexBlock") << "Too many Vertex, fnVertex = " 
                                    << list_->size();
        break;
      }
      vhtm::Vertex vertex;
      vertex.x          = v.x();
      vertex.y          = v.y();
      vertex.z          = v.z();
      vertex.xErr       = v.xError();
      vertex.yErr       = v.yError();
      vertex.zErr       = v.zError();
      vertex.rho        = v.position().rho();
      vertex.chi2       = v.chi2();
      vertex.ndf        = v.ndof();
      vertex.ntracks    = (v.ndof() + 3)/2;
      //vertex.ntracks    = static_cast<int>(v.tracksSize());
      //vertex.ntracksw05 = v.nTracks(0.5); // number of tracks in the vertex with weight above 0.5
      //vertex.sumPt      = v.p4().pt();
      vertex.isfake     = v.isFake();
      vertex.isvalid    = v.isValid();

      if (verbosity_)
	edm::LogInfo("VertexBlock") << std::setw(7) << list_->size()
				    << std::setw(7) << vertex.x 
				    << std::setw(7) << vertex.y 
				    << std::setw(7) << vertex.y
				    << std::setw(7) << vertex.rho
				    << std::setw(9) << vertex.chi2
				    << std::setw(9) << vertex.ndf
				    << std::setw(7) << vertex.ntracks;
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
