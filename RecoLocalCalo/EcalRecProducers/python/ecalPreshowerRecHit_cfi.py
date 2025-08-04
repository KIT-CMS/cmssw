import FWCore.ParameterSet.Config as cms

# Ecal Preshower rechit producer
ecalPreshowerRecHit = cms.EDProducer("ESRecHitProducer",
                                     ESrechitCollection = cms.string('EcalRecHitsES'),
                                     ESdigiCollection = cms.InputTag("ecalPreshowerDigis"),
                                     algo = cms.string("ESRecHitWorker"),
                                     ESRecoAlgo = cms.int32(0)
)
##
## Modify for the tau embedding methods cleaning step
##
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding
from TauAnalysis.MCEmbeddingTools.Cleaning_RECO_cff import common_parameters
tau_embedding.toReplaceWith(ecalPreshowerRecHit, cms.EDProducer("EcalRecHitColCleaner",
    oldCollection = cms.VInputTag(cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES","SELECT")),
    **common_parameters
))