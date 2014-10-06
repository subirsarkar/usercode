#include <memory>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "AnalysisSpace/PattupleMaker/plugins/TauUserEmbedded.h"

TausUserEmbedded::TausUserEmbedded(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
  tauToken_(consumes<pat::TauCollection>(tauTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_))
{
  produces<pat::TauCollection>("");
}

TausUserEmbedded::~TausUserEmbedded()
{
}
void TausUserEmbedded::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::TauCollection> tausHandle;
  iEvent.getByToken(tauToken_, tausHandle);
  const pat::TauCollection* taus = tausHandle.product();

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexToken_, vertexHandle);
  const reco::VertexCollection* vertexes = vertexHandle.product();

  std::auto_ptr<pat::TauCollection> tausUserEmbeddedColl(new pat::TauCollection()) ;
  for (auto it = taus->begin(); it != taus->end(); ++it) {
    pat::Tau pTau((*it));

    float dZPV = vertexes->size() > 0 
      ? std::abs(pTau.vertex().z() - (*vertexes)[0].position().z()) : -99;
    pTau.addUserFloat("dzWrtPV", dZPV);

    tausUserEmbeddedColl->push_back(pTau);
  }
  iEvent.put(tausUserEmbeddedColl);
}

void TausUserEmbedded::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TausUserEmbedded);
