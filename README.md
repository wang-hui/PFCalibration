# PFCalibration  
For Single Pion sample Generation:  

(For tcsh)  
$ setenv SCRAM_ARCH slc6_amd64_gcc530 </br>
(For bash)</br>
$ export SCRAM_ARCH=slc6_amd64_gcc530</br>
$ cmsrel CMSSW_8_1_0_pre16</br>
$ cd CMSSW_8_1_0_pre16/src</br>
$ git clone git@github.com:spandeyehep/PFCalibration.git</br>
$ cd PFCalibration/PFChargedHadronAnalyzer/test/</br>
$ cmsRun PGUnWithGeneration.py
