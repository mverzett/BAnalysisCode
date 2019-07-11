## Instalation:

A clean CMSSW area above 10_2_X. The specific code is being used in 10_5_0. In p
rinciple any area is good.

```
cmsrel CMSSW_10_5_0
cd CMSSW_10_5_0/src
cmsenv
git clone https://github.com/gkaratha/BAnalysisCode
scram b -j 8
```
ready!

Running:
```
cd BAnalysisCode/TriggerAnalyzer/python
cmsRun ConfFile_cfg.py
```

Details of configuration:
all options are in the first 20 lines
IsData= True or False if we run on data/MC
Run="A" takes care of HLT triggers per era A,B,D
Nentries=total evts;  output=output.name; mlog=print a "new event" line per mlog evts;
saveTrk=Saves tracks as collection; NtupleClasses= manages whivh branches are usefull-options: all,auto,class,lite or flat
TrgMuCone=cone to match trg mu with reco (Dr)
UseLowpTe=USE low pT e instead of PF for Kin Fit; LowPtElCollection=Adds low pT e as collection BUT does NOT uses them on Kin Fit; CombinePFLowPtEl=combines PF with low pT e;
LowPtGsfTrkCollection=Adds low pT gsf tracks; PFelCollection=Adds PF e as seperate collection BUT does NOT uses them on kin fit ;
// kinematic/topological cuts for electrons. For low pT seeding BDT values for e1, e2 are tunned from ElXWP. for PF e the code does not read it.
PFLowPtCone is the cone Dr to remove duplicate e when PF low pT is combined
elcuts=dict(El1Pt=1.5,El2Pt=0.5,Dz=0.7,DzeeMax=1.0,El1WP=3,El2WP=-4,PFLowPtCone=0.03)
cuts on trk pT and muon
EtaCut=2.5; TrkPtCut=0.8; MuPtCut=0.5;
RetrieveMuFromTrk=dict(algo=False,maxPtTrk=4)//not tested yet

RecoBtoKLepLep=True; // reconstructs Kll
RecoBtoKstarLepLep=False; //reconstructs K*ll
SkipEventWithNoRecoB=True // skips evts with no reco B
MuonsOnly=False; // if set to true only B in muonic channels will be reconstructed
ElectronsOnly=True; //if set to true only B in electron channels will be reconstruct
ed
addlostTrk=True // adds lost tracks in Packed cand collection
RefitTracksForB="mu"// select if we want refited tracks to be saved. options: mu (only mu), both (mu + e) and none(PF used).Refitted tracks used to calculate B properties.
Bcuts=dict(Prob=0.00000001,Cos=0,MinM=4.5,MaxM=6,MinMll=0,MaxMll=5,PtMin=3.0);
//cuts on reco B
GenRecoMatch=False //when run MC match automatically reco to gen and reconstract only that B
Bdecaymatch=dict(PdgId=521,LepId=11,KId=321,Jtoll=True,DR=0.1)//matching configuration
File=
