import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding


# This modifier sets the correct cuts for the taus decaying into one electron and one muon. More info can be found in TauAnalysis/MCEmbeddingTools.
tau_embedding_emu = cms.ModifierChain(tau_embedding)
