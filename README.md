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
cmsRun PGUnWithGeneration.py
```

