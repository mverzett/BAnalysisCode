import FWCore.ParameterSet.Config as cms
#quick config


IsData=False
Run="D"
Nentries=10000;  output="output_hybrid_test.root"; mlog=1000; 
saveTrk=True; NtupleClasses="auto"; #options: all,auto,class,lite or flat
TrgMuCone=0.03
UseLowpTe=True; LowPtElCollection=False; CombinePFLowPtEl=False;
LowPtGsfTrkCollection=False; PFelCollection=True;
elcuts=dict(El1Pt=0.5,El2Pt=0.5,Dz=1000.7,El1WP=-10.2,El2WP=-10.2,PFLowPtCone=0.03)
TrkPtCut=0.5
RecoBtoKLepLep=True; RecoBtoKstarLepLep=False; SkipEventWithNoRecoB=False
MuonsOnly=False; ElectronsOnly=True; addlostTrk=True
Bcuts=dict(Prob=-0.01,Cos=-10.9,MinM=-4.5,MaxM=60,MinMll=0,MaxMll=5);
GenRecoMatch=True
Bdecaymatch=dict(PdgId=521,LepId=11,KId=321,Jtoll=True,DR=0.1)
File=[
#'/store/data/Run2018B/ParkingBPH5/MINIAOD/PromptReco-v1/000/317/650/00000/321646CB-F76E-E811-91FF-FA163EE936A8.root'
#'/store/mc/RunIIAutumn18MiniAOD/BdToKstar_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13Tev-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v2/80000/F43AD69B-A44E-534B-A195-1C974990B8C6.root'
#'/store/data/Run2018D/ParkingBPH4/MINIAOD/20Mar2019-v1/110000/A336D5D0-93A3-D944-A82A-8FA4165A5B2B.root'
#'/store/data/Run2018D/ParkingBPH1/MINIAOD/20Mar2019-v1/120000/35C993A3-575A-184B-AA17-035629397EDE.root'
#crash #189 evt
#crash #3 evt
#'/store/data/Run2018D/ParkingBPH1/MINIAOD/20Mar2019-v1/120001/44BD5BCD-3EE9-9741-89B6-CB81EE485DA7.root'
'/store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190328_152903/0000/step3_inMINIAODSIM_135.root'
#'file:/afs/cern.ch/work/g/gkaratha/private/SUSYCMG/BtoKppmunu_production/CMSSW_10_2_14/src/BPH-RunIIAutumn18MiniAOD-00001.root'
#'file:/afs/cern.ch/work/g/gkaratha/private/SUSYCMG/HLT/efficiency/Analizer/LostTrackFix/CMSSW_10_5_0/src/step2.root'
#'file:/afs/cern.ch/work/g/gkaratha/private/SUSYCMG/HLT/efficiency/Analizer/LostTrackFix/CMSSW_10_5_0/src/step2.root'
#'/store/mc/RunIIAutumn18MiniAOD/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FE81DFA9-EA47-764C-842D-1F7A31327500.root'
]
############ for debug
 #eventsToProcess=cms.untracked.VEventRange('324293:564:1047993750-324293:564:MAX'),
 #  eventsToProcess=cms.untracked.VEventRange('317696:399:MIN-317696:399:MAX')

##########
electron_container="slimmedElectrons"
if UseLowpTe:
  print "Low pT e collection instead of PF e."
  electron_container="slimmedLowPtElectrons"

electron_container2="slimmedElectrons"
if CombinePFLowPtEl or UseLowpTe or LowPtElCollection or LowPtGsfTrkCollection:
  electron_container2="slimmedLowPtElectrons"

if CombinePFLowPtEl and UseLowpTe:
  print "requested to use low pT as complementary and primary collection -> impossible, using PF as primary"
if CombinePFLowPtEl:
  print "Both collection of e will be used."   

if RecoBtoKLepLep : 
   print "reconstructing B->K(J/psi)ll channel"
if RecoBtoKstarLepLep :
   print "reconstructing B->K*(J/psi)ll->Kpill channel"

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

process.MessageLogger.cerr.FwkReport.reportEvery = mlog #Nentries/100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag, '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(Nentries) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
           File
   ),
   secondaryFileNames=cms.untracked.vstring(
),
 
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
                              electrons    = cms.InputTag(electron_container),
                              lowptElectrons =cms.InputTag(electron_container2),
                              lowptGsftracks=cms.InputTag("lowPtGsfEleGsfTracks"),
                              pfElectrons =cms.InputTag("slimmedElectrons"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              jets = cms.InputTag("slimmedJets"),
                              photons = cms.InputTag("slimmedPhotons"),
                              eleIdMapVeto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                              eleIdMapSoft = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
                              eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
                             eleIdMapTight = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
                              eleIdMapValue = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                              eleBiasedWP = cms.InputTag("lowPtGsfElectronSeedValueMaps","ptbiased","RECO"),
                              eleUnbiasedWP = cms.InputTag("lowPtGsfElectronSeedValueMaps","unbiased","RECO"),
                              #If you want no L1_Seed, write "default" in the first element and the tree will write the value -100
                               Seed=cms.vstring("L1_SingleMu7er1p5","L1_SingleMu8er1p5","L1_SingleMu9er1p5","L1_SingleMu10er1p5","L1_SingleMu12er1p5","L1_SingleMu22"),
                              HLTPath=cms.vstring(n1,n2,n3,n4,n5,n6,n7,n8),
################################NORMALLY USE THIS####################### 
                              triggerresults = cms.InputTag("TriggerResults::HLT"),
                              triggerobjects = cms.InputTag('slimmedPatTrigger',),                   
                              muons=cms.InputTag("slimmedMuons"),
                              met=cms.InputTag("slimmedMETs"),
                              l1seed=cms.InputTag("gtStage2Digis::RECO"),
                              l1met=cms.InputTag('caloStage2Digis','EtSum'), 
                              l1muons=cms.InputTag("gmtStage2Digis","Muon"),
                              l1jets=cms.InputTag('caloStage2Digis','Jet'),                            
                              packed = cms.InputTag("packedGenParticles"),
                              pruned = cms.InputTag("prunedGenParticles"),
                              PFCands=cms.InputTag("packedPFCandidates"),
                              losttracks=cms.InputTag("lostTracks"),

                               RunParameters = cms.PSet(
      Data= cms.bool(IsData),SaveTracks=cms.bool(saveTrk),
      SaveHLT=cms.bool(HLTsave),SaveL1=cms.bool(L1save),
      SaveResultsOnlyIfAPathFired=cms.bool(HLTsave),
      ReconstructBMuMuK=cms.bool(RecoBtoKLepLep),
      ReconstructBMuMuKstar=cms.bool(RecoBtoKstarLepLep),
      MuonPtCutForB=cms.double(1.5),TrackPtCutForB=cms.double(TrkPtCut),
      EtaTrk_Cut=cms.double(2.5),
      #electroncuts
      Electron1PtCut=cms.double(elcuts["El1Pt"]),Electron2PtCut=cms.double(elcuts["El2Pt"]),
      ElectronDzCut=cms.double(elcuts["Dz"]), TrgConeCut=cms.double(TrgMuCone), 
      IsLowpTE=cms.bool(UseLowpTe), MVAEl1Cut=cms.double(elcuts["El1WP"]),
      MVAEl2Cut=cms.double(elcuts["El2WP"]),PointingConstraint=cms.bool(False),
      CosThetaCutPointCons=cms.bool(False), UseClosestVertex=cms.bool(False),
      SkipEventWithNoBToMuMuK=cms.bool(SkipNoKLL),UseBeamspot=cms.bool(False),
      AddeeK=cms.bool(Addel),MLLmax_Cut=cms.double(Bcuts["MaxMll"]),
      MLLmin_Cut=cms.double(Bcuts["MinMll"]),
      MBmin_Cut=cms.double(Bcuts["MinM"]),MBmax_Cut=cms.double(Bcuts["MaxM"]),
      CosThetaCut=cms.double(Bcuts["Cos"]),ProbBMuMuKcut=cms.double(Bcuts["Prob"]), 
      CombineElCol=cms.bool(CombinePFLowPtEl),CombineCone=cms.double(elcuts["PFLowPtCone"]),
      MKstarMin_Cut=cms.double(0.742),MKstarMax_Cut=cms.double(1.042),
      LepTrkExclusionCone=cms.double(0.005),AddLostTracks=cms.bool(addlostTrk),
      RefitTracks=cms.bool(False),RefitMuTracksOnly=cms.bool(True),
      OnlyKee=cms.bool(Onlyel),UsePFeForCos=cms.bool(True),
      SkipEventWithNoBToMuMuKstar=cms.bool(SkipNoKsLL),
      UseDirectlyGenBeeK=cms.bool(GenRecoMatch),DRgenCone=cms.double(Bdecaymatch["DR"]),
      BpdgIdToMatch=cms.int32(Bdecaymatch["PdgId"]),
      LepIdToMatch=cms.int32(Bdecaymatch["LepId"]),
      KIdToMatch=cms.int32(Bdecaymatch["KId"]),
      IsResonantDecayToMatch=cms.bool(Bdecaymatch["Jtoll"]),
      AddLowPtElAsCol=cms.bool(LowPtElCollection),
      AddLowGsfTrkAsCol=cms.bool(LowPtGsfTrkCollection),
      AddPFElAsCol=cms.bool(PFelCollection),
      NtupleOutputClasses=cms.string(NtupleClasses)
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

