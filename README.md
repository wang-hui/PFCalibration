# PFCalibration  
For Single Pion sample Generation:

(For tcsh)  
```
setenv SCRAM_ARCH slc6_amd64_gcc630
```
(For bash)
```
export SCRAM_ARCH=slc6_amd64_gcc630
cmsrel CMSSW_10_0_2
cd CMSSW_10_0_2/src
git clone -b 10_0_2_gen_code git@github.com:spandeyehep/PFCalibration.git
scram b -j 40
cd PFCalibration/PFChargedHadronAnalyzer/test/
cmsenv
cmsRun SinglePi0E10_pythia8_cfi_GEN_SIM.py
```
