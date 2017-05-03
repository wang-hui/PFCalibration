# PFCalibration  
For Single Pion sample Generation:

(For tcsh)  
```
setenv SCRAM_ARCH slc6_amd64_gcc530 </br>
```
(For bash)
```
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_9_0_2
cd CMSSW_9_0_2/src
git clone -b 902_generation_and_calibration git@github.com:spandeyehep/PFCalibration.git
scram b -j 40
cd PFCalibration/PFChargedHadronAnalyzer/test/
cmsenv
cmsRun QCDForPF_14TeV_TuneCUETP8M1_cfi_GEN_SIM.py
```
