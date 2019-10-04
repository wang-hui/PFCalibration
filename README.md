# PFCalibration  
For Single Pion sample Generation:

(For bash)
```
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
git clone -b UltraLegacy_2016 https://github.com/bkansal/PFCalibration.git
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

