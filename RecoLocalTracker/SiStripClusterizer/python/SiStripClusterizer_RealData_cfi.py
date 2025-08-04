import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiStripClusterizer.DefaultClusterizer_cff import *

from RecoLocalTracker.SiStripClusterizer.siStripClusterizer_cfi import siStripClusterizer as _siStripClusterizer
siStripClusters = _siStripClusterizer.clone()

from Configuration.ProcessModifiers.approxSiStripClusters_cff import approxSiStripClusters
from RecoLocalTracker.SiStripClusterizer.SiStripApprox2Clusters_cfi import SiStripApprox2Clusters
SiStripApprox2Clusters.inputApproxClusters = 'SiStripClusters2ApproxClusters'
approxSiStripClusters.toModify(SiStripApprox2Clusters, inputApproxClusters = 'hltSiStripClusters2ApproxClusters')
approxSiStripClusters.toReplaceWith(siStripClusters,SiStripApprox2Clusters)
##
## Modify for the tau embedding methods cleaning step
##
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding
from TauAnalysis.MCEmbeddingTools.Cleaning_RECO_cff import common_parameters
tau_embedding.toReplaceWith(siStripClusters, cms.EDProducer("StripColCleaner",
    oldCollection = cms.VInputTag(cms.InputTag("siStripClusters","","SELECT")),
    **common_parameters
))