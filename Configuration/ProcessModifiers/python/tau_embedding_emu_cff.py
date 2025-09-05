import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding


# This modifier sets the correct cuts for the taus decaying into one electron and one muon. More info can be found in TauAnalysis/MCEmbeddingTools.
_tau_embedding_emu = cms.Modifier()

# As we want to apply the general tau_embedding modifier together with this specific modifier, we create a ModifierChain which can be used in the cmsDriver command
tau_embedding_emu = cms.ModifierChain(tau_embedding, _tau_embedding_emu)
