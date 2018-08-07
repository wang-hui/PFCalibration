# PFCalibration  
For Single Pion sample Generation:

(For tcsh)  
```
setenv SCRAM_ARCH slc6_amd64_gcc530 </br>
```
(For bash)
```
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_1_0_pre16
cd CMSSW_8_1_0_pre16/src
git clone -b PFCalib_eta_dep_check git@github.com:bkansal/PFCalibration.git
scram b -j 40
cd PFCalibration/PFChargedHadronAnalyzer/test/
cmsenv
cmsRun PGUnWithGeneration.py
```
