# PFCalibration  
For Single Pion sample Generation:

(For bash)
```
cmsrel CMSSW_11_0_0_patch1
cd CMSSW_11_0_0_patch1/src
git clone -b for_Run3 https://github.com/bkansal/PFCalibration.git
scram b -j 40
cd PFCalibration/PFChargedHadronAnalyzer/test/
cmsenv
```

Using centerally generated Single pion reco sample for Ultralegacy 2016:
```  
cmsRun myEDAna.py 

(for Crab job submission)
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit crab_step_analyser.py
```
For PFcalibration use PFCalibration/PFChargedHadronAnalyzer/test/calibration_code/calibChris.C  
Note: Please run calibChris.C code on ROOT version 5.xx (it might crash in ROOT version 6.xx)

