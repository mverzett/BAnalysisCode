import FWCore.ParameterSet.Config as cms

unpackedTracksAndVertices = cms.EDProducer(
    'PATTrackAndVertexUnpacker',
    slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    additionalTracks= cms.InputTag("lostTracks"),
    packedCandidates = cms.InputTag("packedPFCandidates")
)

selectedTracks = cms.EDFilter(
    "TrackSelector",
    src = cms.InputTag("unpackedTracksAndVertices"),
    cut = cms.string('')
)

tracks = cms.Sequence(
    unpackedTracksAndVertices *
    selectedTracks
)
