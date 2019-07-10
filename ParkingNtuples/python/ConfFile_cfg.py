import FWCore.ParameterSet.Config as cms
#quick config

selections = cms.PSet(
  electrons = cms.PSet(
    low_pt_collection = cms.InputTag('slimmedLowPtElectrons'),
    low_pt_selection = cms.string('gsfTrack().ptMode() >= 0.5 && passConversionVeto() && electronID("unbiased_seed") >= -4'), #FIXME: add abs(eta) < 2.5
    pf_collection = cms.InputTag('slimmedElectrons'),
    pf_selection = cms.string("passConversionVeto()"), #FIXME: add abs(eta) < 2.5
    cross_cleaning_cone = cms.double(0.03),
    cross_cleaning_dz = cms.double(0.7),
  )
)

IsData=True
Run="B"
Nentries=100;  
output="output_flat.root"; 
mlog=1; 
saveTrk=False; 
NtupleClasses="flat"; #options: all,auto,class,lite or flat

MuonsOnly     = True; 
ElectronsOnly = True; 


TrgMuCone=0.1
UseLowpTe=False; 
LowPtElCollection=False; 
CombinePFLowPtEl=True;
LowPtGsfTrkCollection=False; 
PFelCollection=False;
elcuts=dict(El1Pt=1.5,El2Pt=0.5,Dz=0.7,DzeeMax=1.0,El1WP=3,El2WP=-4,PFLowPtCone=0.03)
EtaCut=2.5; 
TrkPtCut=0.8; 
MuPtCut=0.5;
RetrieveMuFromTrk=dict(algo=False,maxPtTrk=4)
RecoBtoKLepLep=True; 
RecoBtoKstarLepLep=False; 
SkipEventWithNoRecoB=True
addlostTrk=True
RefitTracksForB="mu" #options: mu (only mu), both (mu + e) and none(PF used).Refitted tracks used to calculate B properties.
Bcuts=dict(Prob=0.00000001,Cos=0,MinM=4.5,MaxM=6,MinMll=0,MaxMll=5,PtMin=3.0);
GenRecoMatch=False
Bdecaymatch=dict(PdgId=521,LepId=11,KId=321,Jtoll=True,DR=0.1)
File=[
# '/store/data/Run2018A/ParkingBPH6/MINIAOD/05May2019-v1/260000/6477D465-4909-E34B-A6CE-D7497999E12B.root'
#'/store/user/bainbrid/lowpteleid/BuToKJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_lowpteleid/190328_152903/0000/step3_inMINIAODSIM_130.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2018A/ParkingBPH2/MINIAOD/05May2019-v1/250001/FEECE314-65F3-034B-A6DB-792916EE4EF5.root'
'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'
]
############ for debug evt 1242 run 1 ls 13
 #eventsToProcess=cms.untracked.VEventRange('1:1242:13-1:1242:13'),
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

Addel  = True
Onlyel = False
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
   globaltag='102X_dataRun2_Sep2018Rereco_v1'
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
                  
          ),

)
'''process.selectedPFCandidatesHP = cms.EDFilter("PATPackedCandidateSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("pt() > 0.8 && abs(eta()) < 2.5 && trackHighPurity() > 0")
)'''
#taskB0.add(process.selectedPFCandidatesHP)


## NO NEED!
## from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## dataFormat = DataFormat.MiniAOD
## switchOnVIDElectronIdProducer(process, dataFormat)
## my_id_modules = [
##         'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff', 
##  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
## ]
## for idmod in my_id_modules:
##     setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplesSeq = cms.Sequence(
  #process.egmGsfElectronIDSequence
)

process.load('BAnalysisCode/ParkingNtuples/electrons_cfi')
process.lowptElectronsWithSeed.src = selections.electrons.low_pt_collection
process.lowptElectronsForAnalysis.cut = selections.electrons.low_pt_selection
process.pfElectronsForAnalysis.src = selections.electrons.pf_collection
process.pfElectronsForAnalysis.cut = selections.electrons.pf_selection
process.electronsForAnalysis.drForCleaning = selections.electrons.cross_cleaning_cone
process.electronsForAnalysis.dzForCleaning = selections.electrons.cross_cleaning_dz


process.ntuplesSeq *= process.electrons

process.demo = cms.EDAnalyzer('ParkingNtupleMaker',
                              beamSpot       = cms.InputTag('offlineBeamSpot'),
                              electrons      = cms.InputTag('electronsForAnalysis'),
                              vertices       = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              jets           = cms.InputTag("slimmedJets"),
                              photons        = cms.InputTag("slimmedPhotons"),
                              #If you want no L1_Seed, write "default" in the first element and the tree will write the value -100
                              Seed           = cms.vstring("L1_SingleMu7er1p5","L1_SingleMu8er1p5","L1_SingleMu9er1p5","L1_SingleMu10er1p5","L1_SingleMu12er1p5","L1_SingleMu22"),
                              HLTPath        = cms.vstring(n1,n2,n3,n4,n5,n6,n7,n8),
################################NORMALLY USE THIS####################### 
                              triggerresults = cms.InputTag("TriggerResults::HLT"),
                              triggerobjects = cms.InputTag('slimmedPatTrigger',),                   
                              muons         = cms.InputTag("slimmedMuons"),
                              met           = cms.InputTag("slimmedMETs"),
                              l1seed        = cms.InputTag("gtStage2Digis::RECO"),
                              l1met         = cms.InputTag('caloStage2Digis','EtSum'), 
                              l1muons       = cms.InputTag("gmtStage2Digis","Muon"),
                              l1jets        = cms.InputTag('caloStage2Digis','Jet'),                            
                              packed        = cms.InputTag("packedGenParticles"),
                              pruned        = cms.InputTag("prunedGenParticles"),
                              PFCands       = cms.InputTag("packedPFCandidates"),
                              losttracks    = cms.InputTag("lostTracks"),

                               RunParameters = cms.PSet(
      Data= cms.bool(IsData),
      SaveTracks=cms.bool(saveTrk),
      SaveHLT=cms.bool(HLTsave),
      SaveL1=cms.bool(L1save),
      SaveResultsOnlyIfAPathFired=cms.bool(HLTsave),
      ReconstructBMuMuK=cms.bool(RecoBtoKLepLep),
      ReconstructBMuMuKstar=cms.bool(RecoBtoKstarLepLep),
      MuonMinPtCut = cms.double(1.),
      MuonMaxPtCut = cms.double(1.),
      MuonPtCutForB=cms.double(MuPtCut),
      RetrieveMuFromTrk=cms.bool(RetrieveMuFromTrk["algo"]),
      maxPtTrk=cms.double(RetrieveMuFromTrk["maxPtTrk"]),
      TrackPtCutForB=cms.double(TrkPtCut),
      EtaTrk_Cut=cms.double(EtaCut),
      #electroncuts
      Electron1PtCut=cms.double(elcuts["El1Pt"]),
      Electron2PtCut=cms.double(elcuts["El2Pt"]),
      ElectronDzCut=cms.double(elcuts["Dz"]), 
      DzeeMaxCut=cms.double(elcuts["DzeeMax"]),
      TrgConeCut=cms.double(TrgMuCone), 
      IsLowpTE=cms.bool(UseLowpTe), 
      MVAEl1Cut=cms.double(elcuts["El1WP"]),
      MVAEl2Cut=cms.double(elcuts["El2WP"]),
      PointingConstraint=cms.bool(False),
      CosThetaCutPointCons=cms.bool(False), 
      UseClosestVertex=cms.bool(False),
      SkipEventWithNoBToMuMuK=cms.bool(SkipNoKLL),
      UseBeamspot=cms.bool(False),
      UsePFeForCos=cms.bool(True),
      AddeeK=cms.bool(Addel),
      MLLmax_Cut=cms.double(Bcuts["MaxMll"]),
      MLLmin_Cut=cms.double(Bcuts["MinMll"]),
      MBmin_Cut=cms.double(Bcuts["MinM"]),
      MBmax_Cut=cms.double(Bcuts["MaxM"]),
      PtBminCut=cms.double(Bcuts["PtMin"]),
      CosThetaCut=cms.double(Bcuts["Cos"]),
      ProbBMuMuKcut=cms.double(Bcuts["Prob"]), 
      CombineElCol=cms.bool(CombinePFLowPtEl),
      CombineCone=cms.double(elcuts["PFLowPtCone"]),
      MKstarMin_Cut=cms.double(0.742),
      MKstarMax_Cut=cms.double(1.042),
      LepTrkExclusionCone=cms.double(0.005),
      AddLostTracks=cms.bool(addlostTrk),
      RefitTracks=cms.string(RefitTracksForB),
      OnlyKee=cms.bool(Onlyel),
      SkipEventWithNoBToMuMuKstar=cms.bool(SkipNoKsLL),
      UseDirectlyGenBeeK=cms.bool(GenRecoMatch),
      DRgenCone=cms.double(Bdecaymatch["DR"]),
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
process.ntuplesSeq *= process.demo

process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
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
   #process.egmGsfElectronIDSequence   
   process.ntuplesSeq
   #+process.hlTrigReport
  
   )
   
#process.endjob=cms.EndPath(process.fevt)

