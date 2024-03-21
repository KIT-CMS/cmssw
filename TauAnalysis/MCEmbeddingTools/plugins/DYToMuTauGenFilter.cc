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

class DYToMuTauGenFilter : public edm::stream::EDFilter<> {
public:
  explicit DYToMuTauGenFilter(const edm::ParameterSet &);
  ~DYToMuTauGenFilter() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  bool filter(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  bool leptondecay(const reco::Candidate *d, int depth);
  bool muondecay(const reco::Candidate *d, int depth);

  edm::InputTag inputTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollection_;
  edm::Handle<reco::GenParticleCollection> gen_handle;
};

DYToMuTauGenFilter::DYToMuTauGenFilter(const edm::ParameterSet &iConfig) {
  inputTag_ = iConfig.getParameter<edm::InputTag>("inputTag");
  genParticleCollection_ = consumes<reco::GenParticleCollection>(inputTag_);
}

DYToMuTauGenFilter::~DYToMuTauGenFilter() {}

bool DYToMuTauGenFilter::leptondecay(const reco::Candidate *d, int depth) {
  // returns true if the event has an leptonic decay
  bool check = false;
  if (d->status() != 1) {
    if (d->numberOfDaughters() == 3) {
      if (std::abs(d->daughter(0)->pdgId()) == 14 || std::abs(d->daughter(1)->pdgId()) == 14 ||
          std::abs(d->daughter(2)->pdgId()) == 14 || std::abs(d->daughter(0)->pdgId()) == 12 ||
          std::abs(d->daughter(1)->pdgId()) == 12 || std::abs(d->daughter(2)->pdgId()) == 12)
        check = true;
    } else if (d->numberOfDaughters() > 3)
      return false;
    if (d->numberOfDaughters() < 4) {
      for (unsigned int k = 0; k < d->numberOfDaughters(); k++) {
        int new_depth = depth + 1;
        if (leptondecay(d->daughter(k), new_depth) == true)
          check = true;
      }
    }
  }
  return check;
}
bool DYToMuTauGenFilter::muondecay(const reco::Candidate *d, int depth) {
  // returns true if the event has an muon decay
  bool check = false;
  if (d->status() != 1) {
    if (d->numberOfDaughters() == 3) {
      if (std::abs(d->daughter(0)->pdgId()) == 14 || std::abs(d->daughter(1)->pdgId()) == 14 ||
          std::abs(d->daughter(2)->pdgId()) == 14)
        check = true;
    } else if (d->numberOfDaughters() > 3)
      return false;
    if (d->numberOfDaughters() < 4) {
      for (unsigned int k = 0; k < d->numberOfDaughters(); k++) {
        int new_depth = depth + 1;
        if (muondecay(d->daughter(k), new_depth) == true)
          check = true;
      }
    }
  }
  return check;
}

bool DYToMuTauGenFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  iEvent.getByToken(genParticleCollection_, gen_handle);

  for (unsigned int i = 0; i < gen_handle->size(); i++) {
    const reco::GenParticle gen_particle = (*gen_handle)[i];
    // Check if Z Boson decayed into two leptons
    if (gen_particle.pdgId() == 23 && gen_particle.numberOfDaughters() == 2) {
      // Check if daugther particles are taus
      // From Generator: Mu.Pt > 18 && Had.Pt > 25 && Mu.Eta < 2.1
      if (std::abs(gen_particle.daughter(0)->pdgId()) == 15 && std::abs(gen_particle.daughter(0)->eta()) < 2.1 &&
          gen_particle.daughter(0)->pt() > 25 && gen_particle.daughter(1)->pt() > 18) {
        bool had_1 = !DYToMuTauGenFilter::leptondecay(gen_particle.daughter(0), 1);
        bool mu_2 = DYToMuTauGenFilter::muondecay(gen_particle.daughter(1), 1);
        bool had_2 = !DYToMuTauGenFilter::leptondecay(gen_particle.daughter(1), 1);
        bool mu_1 = DYToMuTauGenFilter::muondecay(gen_particle.daughter(0), 1);

        std::cout << had_1 << " & " << mu_2 << " / " << had_2 << " & " << mu_1 << " |" << std::endl;
        if ((had_1 && mu_2) || (had_2 && mu_1)) {
          std::cout << "Hadronic Decay Check was positive" << std::endl;
          return true;
        }
        std::cout << "Hadronic Decay Check was negative" << std::endl;
      }
      return false;
    }
  }
  return false;
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void DYToMuTauGenFilter::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void DYToMuTauGenFilter::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DYToMuTauGenFilter::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  //  Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DYToMuTauGenFilter);
