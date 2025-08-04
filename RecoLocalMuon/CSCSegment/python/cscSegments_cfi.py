import FWCore.ParameterSet.Config as cms

from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmSK_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmTC_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmDF_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmST_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmRU_cfi import *

from RecoLocalMuon.CSCSegment.cscSegmentProducer_cfi import cscSegmentProducer as _cscSegmentProducer
cscSegments = _cscSegmentProducer.clone(
    # Define input
    inputObjects = cms.InputTag("csc2DRecHits"),
    # Choice of the building algo: 1 SK, 2 TC, 3 DF, 4 ST, 5 RU, ...
    algo_type = cms.int32(5),
    # std::vector<edm::ParameterSet>
    algo_psets = cms.VPSet(
        cms.PSet(
            CSCSegAlgoSK
        ), 
        cms.PSet(
            CSCSegAlgoTC
        ), 
        cms.PSet(
            CSCSegAlgoDF
        ), 
        cms.PSet(
            CSCSegAlgoST
        ),
        cms.PSet(
            CSCSegAlgoRU
        )

     )
)


##
## Modify for the tau embedding methods cleaning step
##
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding
from TauAnalysis.MCEmbeddingTools.Cleaning_RECO_cff import common_parameters
tau_embedding.toReplaceWith(cscSegments, cms.EDProducer("CSCSegmentColCleaner",
    oldCollection = cms.VInputTag(cms.InputTag("cscSegments","","SELECT")),
    **common_parameters
))