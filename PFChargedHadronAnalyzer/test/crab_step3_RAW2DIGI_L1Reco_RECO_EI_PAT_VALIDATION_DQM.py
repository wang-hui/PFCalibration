from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'PGun_step3_RECO_1002_2_200_Feb_13'
config.General.workArea = 'crab_projects'

#optional
#config.General.transferOutputs
#config.General.transferLogs
#config.General.failureLimit = 

#Expert use
#config.General.instance
#config.General.activity

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_EI_PAT_VALIDATION_DQM.py'
config.JobType.psetName = 'step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py'
config.JobType.outputFiles = ['step3.root']
#config.JobType.eventsPerLumi = 2000

config.section_("Data")
#config.Data.inputDataset = '/Single_Pion_gun_13TeV_pythia8/Fall14DR73-NoPU_MCRUN2_73_V9-v1/GEN-SIM-RAW-RECO'
#config.Data.primaryDataset = ''
#config.Data.splitting = 'EventBased'
config.Data.userInputFiles = open('/afs/cern.ch/user/s/spandey/work/public/PF_cal/10_0_2/CMSSW_10_0_2/src/PFCalibration/PFChargedHadronAnalyzer/test/step2_file_list.txt').readlines()
config.Data.ignoreLocality = True
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 5000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = '' default for the moment
#config.Data.outLFN = '/home/spandey/t3store/PF_PGun'
config.Data.outLFNDirBase = '/store/user/spandey/step3/PGun_step3_RECO_1002_2_200_Feb_13/'

config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.blacklist = ['T3_US_UCR', 'T3_US_UMiss']
config.Site.whitelist = ['T2_CH_CERN','T2_KR_KNU']

#config.section_("User")
#config.section_("Debug")
