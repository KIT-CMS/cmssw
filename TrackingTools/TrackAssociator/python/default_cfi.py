import FWCore.ParameterSet.Config as cms

# Default Parameters
#
#   Purpose: extraction of energy deposition and muon matching information
#            for a minimum ionizing particle. Up to 5x5 energy for ECAL 
#            and HCAL should be available.
# 
TrackAssociatorParameterBlock = cms.PSet(
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        useGEM = cms.bool(False),
        GEMSegmentCollectionLabel = cms.InputTag("gemSegments"),
        useME0 = cms.bool(False),
        ME0SegmentCollectionLabel = cms.InputTag("me0Segments"),
        dRHcal = cms.double(9999.0),
        dREcal = cms.double(9999.0),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
	trajectoryUncertaintyTolerance = cms.double(-1.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        accountForTrajectoryChangeCalo = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        preselectMuonTracks = cms.bool(False),
        useCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
	usePreshower = cms.bool(False),
	dRPreshowerPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        RPCHitCollectionLabel = cms.InputTag("rpcRecHits"),
        GEMHitCollectionLabel = cms.InputTag("gemRecHits"),
        ME0HitCollectionLabel = cms.InputTag("me0RecHits"),
        useHcal = cms.bool(True)
    )
)
TrackAssociatorParameters = cms.PSet(
    muonMaxDistanceSigmaX = cms.double(0.0),
    muonMaxDistanceSigmaY = cms.double(0.0),
    CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
    useGEM = cms.bool(False),
    GEMSegmentCollectionLabel = cms.InputTag("gemSegments"),
    useME0 = cms.bool(False),
    ME0SegmentCollectionLabel = cms.InputTag("me0Segments"),
    dRHcal = cms.double(9999.0),
    dREcal = cms.double(9999.0),
    CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
    useEcal = cms.bool(True),
    dREcalPreselection = cms.double(0.05),
    HORecHitCollectionLabel = cms.InputTag("horeco"),
    dRMuon = cms.double(9999.0),
    trajectoryUncertaintyTolerance = cms.double(-1.0),
    propagateAllDirections = cms.bool(True),
    muonMaxDistanceX = cms.double(5.0),
    muonMaxDistanceY = cms.double(5.0),
    useHO = cms.bool(True),
    accountForTrajectoryChangeCalo = cms.bool(False),
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
    EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    dRMuonPreselection = cms.double(0.2),
    usePreshower = cms.bool(False),
    dRHcalPreselection = cms.double(0.2),
    useMuon = cms.bool(True),
    preselectMuonTracks = cms.bool(False),
    useCalo = cms.bool(False),
    EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    truthMatch = cms.bool(False),
    HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
    RPCHitCollectionLabel = cms.InputTag("rpcRecHits"),
    GEMHitCollectionLabel = cms.InputTag("gemRecHits"),
    ME0HitCollectionLabel = cms.InputTag("me0RecHits"),
    useHcal = cms.bool(True)
)
##
## Modify for the tau embedding methods cleaning step
##
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding
tau_embedding.toModify(TrackAssociatorParameterBlock.TrackAssociatorParameters,
    CSCSegmentCollectionLabel = cms.InputTag("cscSegments", "", "SELECT"),
    CaloTowerCollectionLabel = cms.InputTag("towerMaker", "", "SELECT"),
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments", "", "SELECT"),
    EBRecHitCollectionLabel = cms.InputTag("ecalRecHit", "EcalRecHitsEB", "SELECT"),
    EERecHitCollectionLabel = cms.InputTag("ecalRecHit", "EcalRecHitsEE", "SELECT"),
    HBHERecHitCollectionLabel = cms.InputTag("hbhereco", "", "SELECT"),
    HORecHitCollectionLabel = cms.InputTag("horeco", "", "SELECT"),
    ME0HitCollectionLabel = cms.InputTag("me0RecHits", "", "SELECT"),
    ME0SegmentCollectionLabel = cms.InputTag("me0Segments", "", "SELECT"),
    RPCHitCollectionLabel = cms.InputTag("rpcRecHits", "", "SELECT"),
    usePreshower = cms.bool(True)
)