#include <algorithm>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"
#include "TClonesArray.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/GenEventBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

GenEventBlock::GenEventBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genEventTag_(iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfoInputTag", edm::InputTag("generator"))),
  storePDFWeights_(iConfig.getUntrackedParameter<bool>("StorePDFWeights", true)),
  pdfWeightsTag_(iConfig.getUntrackedParameter<edm::InputTag>("PDFWeightsInputTag", edm::InputTag("pdfWeights","cteq66"))),
  genEventToken_(consumes<GenEventInfoProduct>(genEventTag_)),
  pdfWeightsToken_(consumes< std::vector<double> >(pdfWeightsTag_))
{
}
GenEventBlock::~GenEventBlock() {
  delete list_;
  delete pdfWeights_;
}
void GenEventBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  list_ = new std::vector<vhtm::GenEvent>();
  tree->Branch("GenEvent", "std::vector<vhtm::GenEvent>", &list_, 32000, 2);

  pdfWeights_ = new std::vector<double>();
  tree->Branch("pdfWeights", "vector<double>", &pdfWeights_);
}
void GenEventBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector
  list_->clear();

  // Clear the independent vector
  pdfWeights_->clear();

  if (!iEvent.isRealData()) {
    // Create Event Object
    vhtm::GenEvent genEvent;

    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByToken(genEventToken_, genEvtInfoProduct);

    if (genEvtInfoProduct.isValid()) {
      edm::LogInfo("GenEventBlock") << "Success. Obtained GenEventInfoProduct for label: "
                                    << genEventTag_;
      genEvent.processID = genEvtInfoProduct->signalProcessID();
      genEvent.ptHat     = genEvtInfoProduct->hasBinningValues()
                         ? genEvtInfoProduct->binningValues()[0] : 0.;
    }
    else {
      edm::LogError("GenEventBlock") << "Error! Failed to get GenEventInfoProduct for label: "
                                     << genEventTag_;
    }
    // PDF Weights Part
    if (storePDFWeights_) {
      edm::Handle<std::vector<double> > pdfWeightsHandle;
      iEvent.getByToken(pdfWeightsToken_, pdfWeightsHandle);

      if (pdfWeightsHandle.isValid()) {
        edm::LogInfo("GenEventBlock") << "Success. Obtained PDF handle for label: "
                                      << pdfWeightsTag_;
        copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), genEvent.pdfWeights.begin());
        copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), pdfWeights_->begin());
      }
      else {
        edm::LogError("GenEventBlock") << "Error! Failed to get PDF handle for label: "
                                       << pdfWeightsTag_;
      }
    }
    list_->push_back(genEvent);
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
