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
from PhysicsTools.NanoAOD.common_cff import ExtVar, Var

# Keep the selected muons as rows in their own table. Their four-momenta and
# other properties are read directly from the persisted pat::Muon objects.
embeddingMuonTable = cms.EDProducer(
    "SimplePATMuonFlatTableProducer",
    src=cms.InputTag("selectedMuonsForEmbedding"),
    cut=cms.string(""),
    name=cms.string("TauEmbeddingMuon"),
    doc=cms.string("Muons selected for embedding"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        pt=Var("pt", float, doc="muon transverse momentum", precision=10),
        eta=Var("eta", float, doc="muon pseudorapidity", precision=10),
        phi=Var("phi", float, doc="muon azimuthal angle", precision=10),
        mass=Var("mass", float, doc="muon mass", precision=10),
        charge=Var("charge", int, doc="muon charge"),
        dB=Var("dB()", float, doc="muon impact parameter"),
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
