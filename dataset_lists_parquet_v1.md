# Samples description

This README file lists the full path names of the input datasets that were used to produce the parquet file v1.

See cross sections in parquet files dictionaries

## 2022 data

### pre EE+ leak

- [/EGamma/Run2022C-22Sep2023-v1/NANOAOD](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2022C-22Sep2023-v1%2FNANOAOD&instance=prod/global)
- [/EGamma/Run2022D-22Sep2023-v1/NANOAOD](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2022D-22Sep2023-v1%2FNANOAOD&instance=prod/global)

## post EE+ leak

- [/EGamma/Run2022E-22Sep2023-v1/NANOAOD](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2022E-22Sep2023-v1%2FNANOAOD&instance=prod/global)
- [/EGamma/Run2022F-22Sep2023-v1/NANOAOD](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2022F-22Sep2023-v1%2FNANOAOD&instance=prod/global)
- [/EGamma/Run2022G-22Sep2023-v2/NANOAOD](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2022G-22Sep2023-v2%2FNANOAOD&instance=prod/global)



## Non resonant HH signal samples

### Samples ready and ready to process

- The `ZHH_HHto2B2G`, `VBFHHto2B2G`, `WHH_HHto2B2G` on the lists bellow

```
dasgoclient -query=dataset=/*HHto2B2G*/Run3Summer22EENanoAODv12*/*
dasgoclient -query=dataset=/*HHto2B2G*/Run3Summer22NanoAODv12*/*
```

- The BSM couplings for `GluGlutoHHto2B2G` seem to be stuck

## Resonant signal samples

- spin 0 
    - masses: [250, 260, 270, 280, 300, 350, 450, 550, 600, 650, 700, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 4000, 5000]   
    - Run3Summer22postEE: TAG = Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3
    - Run3Summer22preEE: TAG = Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3
- spin 2 
    - masses: [250, 260, 270, 280, 300, 350, 450, 550, 600, 650, 700, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 4000, 5000]  
    - Run3Summer22postEE: TAG = Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3 
    - Run3Summer22preEE: TAG = Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3

- MX = 900 seems missing 

```
/GluGlutoBulkGravitontoHHto2B2G_M-$MASS_narrow_TuneCP5_13p6TeV_madgraph-pythia8/$TAG/NANOAODSIM
```

```
/GluGlutoRadiontoHHto2B2G_M-1000_narrow_TuneCP5_13p6TeV_madgraph-pythia8/$TAG/NANOAODSIM
```

## Background samples

### Non-resonant Background

- GammaGamma+Jets SHERPA (GGJets)
  - preEE: [/GG-Box-3Jets_MGG-80_13p6TeV_sherpa/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGG-Box-3Jets_MGG-80_13p6TeV_sherpa%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2%2FNANOAODSIM&instance=prod/global)
  - postEE: [/GG-Box-3Jets_MGG-80_13p6TeV_sherpa/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGG-Box-3Jets_MGG-80_13p6TeV_sherpa%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global)

- Photon+Jet pt=20-40 GeV PYTHIA (GJetPt20To40)
  - preEE: [/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2%2FNANOAODSIM&instance=prod/global)
  - postEE: [/GJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGJet_PT-20to40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global)

- Photon+Jet pt>40 GeV PYTHIA (GJetPt40)
  - preEE: [/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2%2FNANOAODSIM&instance=prod/global)
  - postEE: [/GJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGJet_PT-40_DoubleEMEnriched_MGG-80_TuneCP5_13p6TeV_pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global)


### Resonant Background

- Gluon fusion production of single H, H->gamma gamma aMC@NLO (GluGluHToGG)
  - preEE: [/GluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2%2FNANOAODSIM&instance=prod/global)
  - postEE: [/GluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FGluGluHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global) 

- VH production of single H, H->gamma gamma aMC@NLO (VHToGG)
  - preEE: [/VHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FVHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2%2FNANOAODSIM&instance=prod/global)
  - postEE: [/VHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v5/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FVHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v5%2FNANOAODSIM&instance=prod/global)

- VBF production of single H, H->gamma gamma aMC@NLO (VBFHToGG)
  - preEE: [/VBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FVBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8%2FRun3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3%2FNANOAODSIM&instance=prod/global)
  - postEE: [/VBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FVBFHtoGG_M-125_TuneCP5_13p6TeV_amcatnlo-pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global)

- ttH production of single H, H->gamma gamma aMC@NLO (ttHToGG)
  - preEE: [/ttHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM]()
  - postEE: [/ttHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM](https://cmsweb.cern.ch/das/request?input=dataset%3D%2FttHtoGG_M-125_TuneCP5_13p6TeV_amcatnloFXFX-madspin-pythia8%2FRun3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2%2FNANOAODSIM&instance=prod/global)


