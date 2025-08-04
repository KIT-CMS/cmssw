import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding

# This modifier is to change tau embedding method to perform mu->e embedding. More info can be found in TauAnalysis/MCEmbeddingTools.
tau_embedding_mu_to_e = cms.ModifierChain(tau_embedding)