# PFCalibration  
For Single Pion sample Generation:

(For tcsh)  
```
setenv SCRAM_ARCH slc6_amd64_gcc530 </br>
```
(For bash)
```
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_9_1_0_3
cd CMSSW_9_1_0_pre3/src
git clone -b For_91X_devel git@github.com:spandeyehep/PFCalibration.git
scram b -j 40
cmsenv
```
