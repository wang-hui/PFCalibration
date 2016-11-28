# PFCalibration
For Single Pion sample Generation:

#For tcsh
$ setenv SCRAM_ARCH slc6_amd64_gcc530
#For bash
$ export SCRAM_ARCH=slc6_amd64_gcc530
$ cmsrel CMSSW_8_1_0_pre16
$ cd CMSSW_8_1_0_pre16/src
$ git clone git@github.com:spandeyehep/PFCalibration.git
$ cd PFCalibration/PFChargedHadronAnalyzer/test/
$ cmsRun PGUnWithGeneration.py
