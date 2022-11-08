/** \class TrackerCleaner
 *
 *
 * \author Stefan Wayand;
 *         Christian Veelken, LLR
 *
 *
 *
 *
 *
 */

#ifndef TauAnalysis_MCEmbeddingTools_TrackerCleaner_H
#define TauAnalysis_MCEmbeddingTools_TrackerCleaner_H

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"


#include <iostream>
#include <map>
#include <string>

template <typename T> class TrackerCleaner : public edm::stream::EDProducer<>
{
public:
    explicit TrackerCleaner(const edm::ParameterSet &);
    ~TrackerCleaner() override;

private:
    void produce(edm::Event &, const edm::EventSetup &) override;

    const edm::EDGetTokenT<edm::View<pat::Muon>> mu_input_;
    typedef edmNew::DetSetVector<T> TrackClusterCollection;

    std::map<std::string, edm::EDGetTokenT<TrackClusterCollection>> inputs_;

    bool match_rechit_type(const TrackingRecHit &murechit);
};

template <typename T> TrackerCleaner<T>::TrackerCleaner(const edm::ParameterSet &iConfig) : mu_input_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonCollection")))

{
    std::vector<edm::InputTag> inCollections = iConfig.getParameter<std::vector<edm::InputTag>>("oldCollection");
    for (auto inCollection : inCollections)
    {
        inputs_[inCollection.instance()] = consumes<TrackClusterCollection>(inCollection);
        produces<TrackClusterCollection>(inCollection.instance());
    }
}

template <typename T> TrackerCleaner<T>::~TrackerCleaner()
{
    // nothing to be done yet...
}

template <typename T> void TrackerCleaner<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{

    std::cout << "Start TrackerCleaner" << std::endl;

    using namespace edm;

    edm::Handle<edm::View<pat::Muon>> muonHandle;
    iEvent.getByToken(mu_input_, muonHandle);
    edm::View<pat::Muon> muons = *muonHandle;

    // std::cout << "Considered type: " << typeid(TrackClusterCollection).name() << std::endl;

    for (auto input_ : inputs_)
    {

        edm::Handle<TrackClusterCollection> inputClusters;
        iEvent.getByToken(input_.second, inputClusters);

        std::vector<bool> vetodClusters;

        vetodClusters.resize(inputClusters->dataSize(), false);

        for (edm::View<pat::Muon>::const_iterator iMuon = muons.begin(); iMuon != muons.end(); ++iMuon)
        {
            if (!iMuon->isGlobalMuon())
                continue;
            const reco::Track *mutrack = iMuon->globalTrack().get();
            //  reco::Track *mutrack = new reco::Track(*(iMuon->innerTrack() ));
            for (trackingRecHit_iterator hitIt = mutrack->recHitsBegin(); hitIt != mutrack->recHitsEnd(); ++hitIt)
            {
                const TrackingRecHit &murechit = **hitIt;
                if (!(murechit).isValid())
                    continue;

                if (match_rechit_type(murechit))
                {
                    // std::cout << "\trawId from murechit: " << murechit.rawId() << ", rawId from iterator: " << (*hitIt)->rawId();
                    // std::cout << ", type: " << (*hitIt)->getType() << ", is valid? " << (*hitIt)->isValid();
                    // std::cout << ", local x, y, z: " << (*hitIt)->localPosition();

                    auto &thit = reinterpret_cast<BaseTrackerRecHit const &>(murechit);
                    auto const &cluster = thit.firstClusterRef();

                    // std::cout << ", firstClusterRef: " << cluster.key() << std::endl;

                    // std::cout << "Before local position" << std::endl;
                    // auto localPos = murechit.globalPosition();
                    // std::cout << "Between" << std::endl;
                    // // auto pos_x = globPos.x(); // float
                    // std::cout << "Local Position : " << localPos << std::endl;
                    // std::cout << "After local position" << std::endl;
                    // std::cout << "&(thit.globalPosition()): " << &(thit.globalPosition()) << std::endl;

                    vetodClusters[cluster.key()] = true;
                }
                auto &thit = reinterpret_cast<BaseTrackerRecHit const &>(murechit);
                if (trackerHitRTTI::isMatched(thit))
                {
                    vetodClusters[reinterpret_cast<SiStripMatchedRecHit2D const&>(murechit).stereoClusterRef().key()] = true;
                }
            }
        }
        /// New output collection for the producer to write out.
        std::unique_ptr<TrackClusterCollection> output(new TrackClusterCollection());

        int idx = 0;
        int set_idx = 0;
        // Loop over Cluster Collection to receive ClusterSets
        for (typename TrackClusterCollection::const_iterator clustSet = inputClusters->begin(); clustSet != inputClusters->end(); ++clustSet)
        {
            set_idx++;
            DetId detIdObject = clustSet->detId();
            typename TrackClusterCollection::FastFiller spc(*output, detIdObject);
            // Loop over Cluster Sets to recieve a detector set
            for (typename edmNew::DetSet<T>::const_iterator clustIt = clustSet->begin(); clustIt != clustSet->end(); ++clustIt)
            {
                if (vetodClusters[idx]) {
                    std::string name = typeid(*clustIt).name();
//                    std::cout << "\t\tDetectorSet: " << name << ", DetectorSet size: " << clustIt->size() << ", clusterSet size: " << clustSet->size()  << ", Detector object: " << detIdObject.rawId() << ", index: " << idx << ", set index: " << set_idx;
                    // if(name.find("Strip") != std::string::npos)
                    // {
                    //    std::cout << ", barycenter: " << (reinterpret_cast<const SiStripCluster*>(clustIt))->barycenter() << std::endl;
                    // }
                    // else if (name.find("Pixel") != std::string::npos)
                    // {
                    //    std::cout << ", x: " << (reinterpret_cast<const SiPixelCluster*>(clustIt))->x() << ", y: " << (reinterpret_cast<const SiPixelCluster*>(clustIt))->y() << std::endl;
                    // }
                }
                else {
                    // if (!vetodClusters[idx-1]) continue; for inverted selction
                    spc.push_back(*clustIt);
                }
                idx++;
            }
        }
        iEvent.put(std::move(output), input_.first);
    }
}
#endif
