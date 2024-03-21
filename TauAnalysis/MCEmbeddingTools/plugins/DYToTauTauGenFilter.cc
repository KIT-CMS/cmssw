#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DYToTauTauGenFilter : public edm::stream::EDFilter<> {
public:
  explicit DYToTauTauGenFilter(const edm::ParameterSet &);
  ~DYToTauTauGenFilter() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  bool filter(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;

  edm::InputTag inputTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollection_;

  edm::Handle<reco::GenParticleCollection> gen_handle;

  // ----------member data ---------------------------
};

DYToTauTauGenFilter::DYToTauTauGenFilter(const edm::ParameterSet &iConfig) {
  inputTag_ = iConfig.getParameter<edm::InputTag>("inputTag");
  genParticleCollection_ = consumes<reco::GenParticleCollection>(inputTag_);
}

DYToTauTauGenFilter::~DYToTauTauGenFilter() {}

bool photoncheck(const reco::Candidate *d, int depth) {
  bool check = false;
  if (d->status() != 1) {
    if (d->numberOfDaughters() == 3) {
      if (std::abs(d->daughter(0)->pdgId()) == 14 || std::abs(d->daughter(1)->pdgId()) == 14 ||
          std::abs(d->daughter(2)->pdgId()) == 14 || std::abs(d->daughter(0)->pdgId()) == 12 ||
          std::abs(d->daughter(1)->pdgId()) == 12 || std::abs(d->daughter(2)->pdgId()) == 12) {
        check = true;
      }
    } else if (d->numberOfDaughters() > 3)
      return false;
    if (d->numberOfDaughters() < 4) {
      for (unsigned int k = 0; k < d->numberOfDaughters(); k++) {
        int new_depth = depth + 1;
        if (photoncheck(d->daughter(k), new_depth) == true)
          check = true;
      }
    }
  }
  return check;
}

bool DYToTauTauGenFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  iEvent.getByToken(genParticleCollection_, gen_handle);

  for (unsigned int i = 0; i < gen_handle->size(); i++) {
    const reco::GenParticle gen_particle = (*gen_handle)[i];
    // Check if Z Boson decayed into two leptons
    if (gen_particle.pdgId() == 23 && gen_particle.numberOfDaughters() == 2) {
      // Check if daugther particles are taus
      if (std::abs(gen_particle.daughter(0)->pdgId()) == 15 && std::abs(gen_particle.daughter(0)->eta()) < 2.3 &&
          std::abs(gen_particle.daughter(1)->eta()) < 2.3 &&
          ((gen_particle.daughter(0)->pt() > 30 && gen_particle.daughter(1)->pt() > 35) ||
           (gen_particle.daughter(0)->pt() > 35 && gen_particle.daughter(1)->pt() > 30))) {
        if (!photoncheck(gen_particle.daughter(0), 1) && !photoncheck(gen_particle.daughter(1), 1)) {
          return true;
        }
      }
      return false;
    }
  }
  return false;
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void DYToTauTauGenFilter::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void DYToTauTauGenFilter::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DYToTauTauGenFilter::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DYToTauTauGenFilter);
