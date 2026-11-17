// -*- C++ -*-
//
// Package:    TauAnalysis/EmbeddingProducer
// Class:      MuMuForEmbeddingSelector
//
/**\class MuMuForEmbeddingSelector MuMuForEmbeddingSelector.cc TauAnalysis/EmbeddingProducer/plugins/MuMuForEmbeddingSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Artur Akhmetshin
//         Created:  Mon, 13 Jun 2016 11:05:32 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//
// class declaration
//

class MuMuForEmbeddingSelector : public edm::stream::EDProducer<> {
public:
  explicit MuMuForEmbeddingSelector(const edm::ParameterSet &);

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate>> ZmumuCandidates_;
  edm::EDGetTokenT<edm::View<pat::MET>> theMETLabel_;
  edm::EDGetTokenT<edm::View<pat::MET>> thePuppiMETLabel_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuMuForEmbeddingSelector::MuMuForEmbeddingSelector(const edm::ParameterSet &iConfig)
    : ZmumuCandidates_(consumes<edm::View<reco::CompositeCandidate>>(
          iConfig.getParameter<edm::InputTag>("ZmumuCandidatesCollection"))) {
  produces<edm::RefVector<pat::MuonCollection>>();
  // Keep references to the selected PAT muons. The source collection must be
  // retained through merging for the NanoAOD table to dereference these refs.
  produces<float>("initialMETEt");
  produces<float>("initialMETphi");
  produces<float>("initialPuppiMETEt");
  produces<float>("initialPuppiMETphi");
  theMETLabel_ = consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("Met"));
  thePuppiMETLabel_ = consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("PuppiMet"));
  // now do what ever other initialization is needed
}

// 
// member functions
//

// ------------ method called to produce the data  ------------
void MuMuForEmbeddingSelector::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace edm;
  edm::Handle<edm::View<reco::CompositeCandidate>> ZmumuCandidatesHandle;
  iEvent.getByToken(ZmumuCandidates_, ZmumuCandidatesHandle);
  edm::View<reco::CompositeCandidate> ZmumuCandidates = *ZmumuCandidatesHandle;
  const reco::CompositeCandidate *chosenZCand = nullptr;
  edm::Handle<edm::View<pat::MET>> met;
  iEvent.getByToken(theMETLabel_, met);
  edm::Handle<edm::View<pat::MET>> puppimet;
  iEvent.getByToken(thePuppiMETLabel_, puppimet);
  for (edm::View<reco::CompositeCandidate>::const_iterator iZCand = ZmumuCandidates.begin();
       iZCand != ZmumuCandidates.end();
       ++iZCand) {
    if (chosenZCand == nullptr) {
      chosenZCand = &(*iZCand);
    } else {
      if (iZCand->mass() > chosenZCand->mass()) {
        chosenZCand = &(*iZCand);
      }
    }
  }

  const auto &daughter0 = chosenZCand->daughter(0)->masterClone().castTo<pat::MuonRef>();
  const auto &daughter1 = chosenZCand->daughter(1)->masterClone().castTo<pat::MuonRef>();

  std::unique_ptr<edm::RefVector<pat::MuonCollection>> prod(new edm::RefVector<pat::MuonCollection>());
  prod->reserve(2);
  prod->push_back(daughter0);
  prod->push_back(daughter1);
  iEvent.put(std::move(prod));

  iEvent.put(std::make_unique<float>(met->at(0).et()), "initialMETEt");
  iEvent.put(std::make_unique<float>(met->at(0).phi()), "initialMETphi");
  iEvent.put(std::make_unique<float>(puppimet->at(0).et()), "initialPuppiMETEt");
  iEvent.put(std::make_unique<float>(puppimet->at(0).phi()), "initialPuppiMETphi");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuMuForEmbeddingSelector::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(MuMuForEmbeddingSelector);
