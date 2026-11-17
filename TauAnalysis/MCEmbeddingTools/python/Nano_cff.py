"""
This config fragment is used to generate NanoAOD samples from tau embedding MiniAOD samples.
It adds an additional table with information about the initially selected muons and the initial event.

As the tau embedding event is a hybrid event and we want to have information about the measured event and the generated tau decay,
Both, the normal NanoAOD sequence for measured events and the NanoAOD sequence for simulated events, have to run.
This and the embeddingTable_seq is configured in the autoNANO dict in PhysicsTools/NanoAOD/python/autoNANO.py
With this one can simply add `NANO:@TauEmbedding` to the `--step` parameter of the cmsDriver command, to get all the additional tables.
The merging step must be carried out beforehand.
To use this config fragment, a cmsDriver command like the following can be used:
```
cmsDriver.py \
	--step NANO:@TauEmbedding \
	--data \
	--scenario pp \
	--eventcontent TauEmbeddingNANOAOD \
	--datatier NANOAODSIM \
	--era ... \
    --conditions ... \
    --filein ... \
    --fileout ...
```
"""

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars, ExtVar, Var
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer

# Preserve the two selected muons as rows, with their complete four-vectors
# and standard identification flags, instead of duplicating their properties
# as leading/trailing event-level float products.
embeddingMuonTable = simplePATMuonFlatTableProducer.clone(
    src=cms.InputTag("selectedMuonsForEmbedding"),
    name=cms.string("TauEmbeddingMuon"),
    doc=cms.string("Muons selected in the original event for embedding"),
    variables=cms.PSet(
        CandVars,
        dB=Var("dB()", float, doc="transverse impact parameter"),
        vx=Var("vx", float, doc="x coordinate of the muon vertex"),
        vy=Var("vy", float, doc="y coordinate of the muon vertex"),
        vz=Var("vz", float, doc="z coordinate of the muon vertex"),
        tightId=Var("passed('CutBasedIdTight')", bool, doc="tight muon ID"),
        # looseId=Var("passed('CutBasedIdLoose')", bool, doc="loose muon ID"),
        # mediumId=Var("isMediumMuon", bool, doc="medium muon ID"),
        # isPFcand=Var("isPFMuon", bool, doc="muon is a PF muon"),
        # isGlobal=Var("isGlobalMuon", bool, doc="muon is a global muon"),
        # isTracker=Var("isTrackerMuon", bool, doc="muon is a tracker muon"),
    ),
)

# This table producer adds information about the initially selected muons and the initial event to the NanoAOD.
# The information comes from the `selectedMuonsForEmbedding` producer, which is created in the selection step.
embeddingTable = cms.EDProducer(
    "GlobalVariablesTableProducer",
    name=cms.string("TauEmbedding"),
    # doc=cms.string("TauEmbedding"),
    variables=cms.PSet(
        initialMETEt=ExtVar(
            cms.InputTag("selectedMuonsForEmbedding", "initialMETEt"),
            float,
            doc="MET Et of selected event",
        ),
        initialMETphi=ExtVar(
            cms.InputTag("selectedMuonsForEmbedding", "initialMETphi"),
            float,
            doc="MET phi of selected event",
        ),
        initialPuppiMETEt=ExtVar(
            cms.InputTag("selectedMuonsForEmbedding", "initialPuppiMETEt"),
            float,
            doc="PuppiMET Et of selected event",
        ),
        initialPuppiMETphi=ExtVar(
            cms.InputTag("selectedMuonsForEmbedding", "initialPuppiMETphi"),
            float,
            doc="PuppiMET phi of selected event",
        ),
    ),
)
# put all sequences and tasks together to create the embedding NanoAOD sequence
embeddingTableTask = cms.Task(embeddingMuonTable, embeddingTable)
embeddingTable_seq = cms.Sequence(embeddingTableTask)
