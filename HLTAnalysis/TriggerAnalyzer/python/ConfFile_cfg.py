import FWCore.ParameterSet.Config as cms
#quick config
IsData=True
Run="A"
output="output_KstarMuMu.root"
RecoBtoKLepLep=False
RecoBtoKstarLepLep=True
SkipEventWithNoRecoB=False
MuonsOnly=True
ElectronsOnly=False
addlostTrk=True
saveTrk=False
Nentries=1000
File=['/store/data/Run2018B/ParkingBPH5/MINIAOD/PromptReco-v1/000/317/650/00000/321646CB-F76E-E811-91FF-FA163EE936A8.root']
############
if RecoBtoKLepLep : 
   print "reconstructing B->Kll channel"
if RecoBtoKstarLepLep :
   print "reconstructing B->K*ll->Kpill channel"

Addel=True
Onlyel=False
if MuonsOnly and not ElectronsOnly:
  Addel=False 
elif ElectronsOnly and not MuonsOnly:
  Onlyel=True
elif ElectronsOnly and MuonsOnly:
   print "warning both chanels will be kept"


SkipNoKLL=False
if RecoBtoKLepLep and SkipEventWithNoRecoB :
   SkipNoKLL=True
#print SkipNoKLL
SkipNoKsLL=False
if RecoBtoKstarLepLep and SkipEventWithNoRecoB :
   SkipNoKsLL=True

if Run=="A":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu8p5_IP3p5" ; n3 ="HLT_Mu10p5_IP3p5" ; 
   n4="HLT_Mu8_IP3"; n5="empty" ; n6="empty" ; n7="empty" ; n8="empty" 
elif Run=="B":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu9_IP5" ; n3 ="HLT_Mu7_IP4" ; 
   n4="HLT_Mu8_IP3"; n5="HLT_Mu12_IP6" ; n6="empty" ; n7="empty" ; n8="empty"
elif Run=="D":
   n1="HLT_Mu9_IP6_part" ; n2="HLT_Mu9_IP5" ; n3="HLT_Mu7_IP4"; n4="HLT_Mu8_IP3"; n5="HLT_Mu12_IP6" ; n6="HLT_Mu9_IP4" ; n7="HLT_Mu8_IP6"; n8="HLT_Mu8_IP5"
else:
   n1="empty" ; n2=n1 ;n3=n1 ; n4=n1 ; n5=n1 ; n6=n1 ; n7=n1 ; n8=n1
   

globaltag='102X_upgrade2018_realistic_v15' 
L1save=False ; HLTsave=False ; HLTfired=False
if IsData:
   print "We have established we Run on data"
   globaltag='101X_dataRun2_Prompt_v11'
   L1save=True ; HLTsave=True ; HLTfired=True
else:
   print "We have established we Run on MC"
print "Run parameters ",globaltag," save L1 ",L1save," HLT ",HLTsave," save ev. if path fired ",HLTfired

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#tracks from pf
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(Nentries) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    File
#for aod comparison
# '/store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/F637E740-F259-E811-B4BB-782BCB3B1A58.root'
#'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_MINIAODSIM_18_03_21/180319_172427/0001/BToKee_MINIAODSIM_1727.root'
#aod comp 2
#'/store/data/Run2018B/ParkingBPH5/MINIAOD/PromptReco-v1/000/317/696/00000/F216061E-AF70-E811-B98A-FA163EC375AA.root'
#'/store/user/tstreble/BToKmm_Pythia/BToKmm_Pythia_MINIAODSIM_18_03_26/180326_093119/0000/BToKmm_MINIAODSIM_99.root'
#fastsim
#"file:/afs/cern.ch/work/g/gkaratha/private/SUSYCMG/BtoKppmunu_production/CMSSW_9_3_6/src/step3_PAT.root"
#'/store/group/cmst3/user/gkaratha/PAT_FastSim_GenBToKEE_NonD_muFilter_ForProbe_try1/CRAB_UserFiles/crab_PAT_FastSim_GenBToKEE_NonD_muFilter_ForProbe_try1/190129_202905/0000/step3_PAT_2.root'
#'/store/mc/RunIIAutumn18MiniAOD/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FE81DFA9-EA47-764C-842D-1F7A31327500.root'
  ),
   secondaryFileNames=cms.untracked.vstring(
),
 # eventsToProcess=cms.untracked.VEventRange('316187:796:MIN-316187:796:MAX'),
 #  eventsToProcess=cms.untracked.VEventRange('317696:399:MIN-317696:399:MAX')
  #eventsToProcess=cms.untracked.VEventRange('317696:252:329504213-317696:252:329504215'),
   inputCommands=cms.untracked.vstring(
                  'keep *',
                  'drop *_ctppsPixelClusters_*_*',
                  
          )

)
'''process.selectedPFCandidatesHP = cms.EDFilter("PATPackedCandidateSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("pt() > 0.8 && abs(eta()) < 2.5 && trackHighPurity() > 0")
)'''
#taskB0.add(process.selectedPFCandidatesHP)


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff', 
 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.demo = cms.EDAnalyzer('TriggerAnalyzerb',
                              beamSpot = cms.InputTag('offlineBeamSpot'),
                              electrons    = cms.InputTag("slimmedElectrons"),
                              vertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              jets = cms.InputTag("slimmedJets"),
                              photons = cms.InputTag("slimmedPhotons"),
                              eleIdMapVeto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                              eleIdMapSoft = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
                              eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
                             eleIdMapTight = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
                              eleIdMapValue = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                              #If you want no L1_Seed, write "default" in the first element and the tree will write the value -100
                               Seed=cms.vstring("L1_SingleMu7er1p5","L1_SingleMu8er1p5","L1_SingleMu9er1p5","L1_SingleMu10er1p5","L1_SingleMu12er1p5","L1_SingleMu22"),
                              HLTPath=cms.vstring(n1,n2,n3,n4,n5,n6,n7,n8),
################################NORMALLY USE THIS####################### 
                              triggerresults = cms.InputTag("TriggerResults::HLT"),
                              triggerobjects = cms.InputTag('slimmedPatTrigger','','RECO'),                   
                              muons=cms.InputTag("slimmedMuons"),
                              met=cms.InputTag("slimmedMETs"),
                              l1seed=cms.InputTag("gtStage2Digis::RECO"),
                              l1met=cms.InputTag('caloStage2Digis','EtSum','RECO'), 
                              l1muons=cms.InputTag("gmtStage2Digis","Muon","RECO"),
                              l1jets=cms.InputTag('caloStage2Digis','Jet','RECO'),                            
                              packed = cms.InputTag("packedGenParticles"),
                              pruned = cms.InputTag("prunedGenParticles"),
                              PFCands=cms.InputTag("packedPFCandidates"),
                              losttracks=cms.InputTag("lostTracks"),
 #                              tracks=cms.InputTag("unpackedTracksAndVertices"),

                               RunParameters = cms.PSet(
      Data= cms.bool(IsData),SaveTracks=cms.bool(saveTrk),
      SaveHLT=cms.bool(HLTsave),SaveL1=cms.bool(L1save),
      SaveResultsOnlyIfAPathFired=cms.bool(HLTsave),
      ReconstructBMuMuK=cms.bool(RecoBtoKLepLep),
      ReconstructBMuMuKstar=cms.bool(RecoBtoKstarLepLep),
      MuonPtCutForB=cms.double(1.5),TrackPtCutForB=cms.double(1.3),
      PointingConstraint=cms.bool(False),CosThetaCutPointCons=cms.bool(False),
      CosThetaCut=cms.double(-1),UseClosestVertex=cms.bool(False),
      ProbBMuMuKcut=cms.double(-0.001),
      SkipEventWithNoBToMuMuK=cms.bool(SkipNoKLL),UseBeamspot=cms.bool(False),
      AddeeK=cms.bool(Addel),MLLmax_Cut=cms.double(5),MLLmin_Cut=cms.double(0),
      MBmin_Cut=cms.double(4.5),
      MBmax_Cut=cms.double(6),EtaTrk_Cut=cms.double(2.5),
      MKstarMin_Cut=cms.double(0.742),MKstarMax_Cut=cms.double(1.042),
      LepTrkExclusionCone=cms.double(0.005),AddLostTracks=cms.bool(addlostTrk),
      RefitTracks=cms.bool(False),RefitMuTracksOnly=cms.bool(True),
      OnlyKee=cms.bool(Onlyel),UsePFeForCos=cms.bool(True),
      SkipEventWithNoBToMuMuKstar=cms.bool(SkipNoKsLL)
  ),
)

process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.hlTrigReport.HLTriggerResults   = cms.InputTag("TriggerResults", "", "HLT")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(output)
                                   )
process.fevt = cms.OutputModule("PoolOutputModule",
   # SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring(#"drop *",
    ),
    fileName = cms.untracked.string("edm_output.root"))

#process.p = cms.Path(process.egmGsfElectronIDSequence)#* process.demo)
process.p = cms.Path(
   process.egmGsfElectronIDSequence   
  # +process.unpackedTracksAndVertices
#+process.SecondaryVerticesFromHighPurityTracks
 #  +process.selectedPFCandidatesHP
   +process.demo
   #+process.hlTrigReport
  
   )
   
#process.endjob=cms.EndPath(process.fevt)
#samples
#'/store/mc/RunIIAutumn18MiniAOD/BuToK_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/FBEB6F2C-3302-9C4E-9E9D-F253120EC027.root'

