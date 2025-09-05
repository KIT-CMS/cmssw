import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.tau_embedding_cff import tau_embedding

# This modifier is to change tau embedding method to perform mu->mu embedding. More info can be found in TauAnalysis/MCEmbeddingTools.
_tau_embedding_mu_to_mu = cms.Modifier()

# As we want to apply the general tau_embedding modifier together with this specific modifier, we create a ModifierChain which can be used in the cmsDriver command
tau_embedding_mu_to_mu =  cms.ModifierChain(tau_embedding, _tau_embedding_mu_to_mu)