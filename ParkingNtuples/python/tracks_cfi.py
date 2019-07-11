import FWCore.ParameterSet.Config as cms

mergedCandidates = cms.EDProducer(
    'PackedCandidateMerger',
    packedSrc = cms.InputTag("packedPFCandidates"),
    lostSrc = cms.InputTag("lostTracks"),
)

selectedTracks = cms.EDFilter(
    "PackedCandidateSelector",
    src = cms.InputTag("mergedCandidates"),
    cut = cms.string('')
)

tracks = cms.Sequence(
    mergedCandidates *
    selectedTracks
)
