import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def ufloat(expr, precision=-1):
    return Var('userFloat("%s")' % expr, float, precision = precision)

muonTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("FIXME"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("Our muon table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        CandVars,
        ),
    externalVariables = cms.PSet(),
)

electronTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("FIXME"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Electron"),
    doc  = cms.string("Our electron table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        pt  = Var('? (userInt("isLowPt") == 1) ? gsfTrack().ptMode()  : pt' ,  float, precision=-1),
        eta = Var('? (userInt("isLowPt") == 1) ? gsfTrack().etaMode() : eta',  float, precision=12),
        phi = Var('? (userInt("isLowPt") == 1) ? gsfTrack().phiMode() : phi',  float, precision=12),
        charge = Var('charge', int),
        ),
    externalVariables = cms.PSet(),
)

trackTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("FIXME"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Track"),
    doc  = cms.string("Our packedCandidate table"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        CandVars,
        ),
    externalVariables = cms.PSet(),
)

bToKMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("FIXME"),
    cut = cms.string(""),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        CandVars,
        mu1Idx = Var('userInt("l1_idx")', int),
        mu2Idx = Var('userInt("l2_idx")', int),
        kIdx   = Var('userInt("k_idx" )', int),
    ) 
)

bToKMuMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("FIXME"),
    cut = cms.string(""),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        CandVars,
        mu1Idx = Var('userInt("l1_idx")', int),
        mu2Idx = Var('userInt("l2_idx")', int),
        kIdx   = Var('userInt("k_idx" )', int),
    ) 
)

bToKEETable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("FIXME"),
    cut = cms.string(""),
    name = cms.string("BToKEE"),
    doc = cms.string("BToKEE Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        CandVars,
        el1Idx = Var('userInt("l1_idx")', int),
        el2Idx = Var('userInt("l2_idx")', int),
        kIdx   = Var('userInt("k_idx" )', int),
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        mll = ufloat('m_ll'),
        cos2D = ufloat('cos_theta_2D'),
    ) 
)

tables = cms.Sequence(
    muonTable + 
    electronTable +
    trackTable +
    bToKMuMuTable +
    bToKEETable
)
