import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding

# This modifier sets the correct cuts for the taus decaying into one jet and one electron. More info can be found in TauAnalysis/MCEmbeddingTools.
tau_embedding_etauh = cms.ModifierChain(tau_embedding)
