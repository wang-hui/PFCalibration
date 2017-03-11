# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --runUnscheduled --conditions auto:phase1_2017_hcaldev -n 10 --era Run2_2017_HCALdev --eventcontent RECOSIM -s RAW2DIGI,L1Reco,RECO,EI --datatier GEN-SIM-RECO --geometry Extended2017dev -n 100 --nThreads 16
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT3',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("root://se01.indiacms.res.in//store/user/spandey/PGun_step2_DIGI_Feb_8/CRAB_UserFiles/crab_PGun_step2_DIGI_Feb_8/170208_131436/0000/step2_pi1-100_2017_realistic_1.root"),
                            secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:step3_pi1-100_2017_realistic.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.RECOSIMoutput.outputCommands.append( "keep *_g4SimHits*_*_*")
process.RECOSIMoutput.outputCommands.append( "keep *_hbheprereco*_*_*")
process.RECOSIMoutput.outputCommands.append( "keep HBHERecHits*_*_*_*")

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '81X_upgrade2017_realistic_v26', '')

# Save HBHEChannelInfo
#process.hbheprereco.saveInfos = cms.bool(True)

#process.hbheprerecoM3 = process.hbheprereco.clone()
#process.hbheprerecoM3.algorithm.__setattr__('useM2',cms.bool(False))
#process.hbheprerecoM3.algorithm.__setattr__('useM3',cms.bool(True))





process.pfChargedHadronAnalyzer = cms.EDAnalyzer(
    "PFChargedHadronAnalyzer",
    PFCandidates = cms.InputTag("particleFlow"),
    PFSimParticles = cms.InputTag("particleFlowSimParticle"),
    EcalPFClusters = cms.InputTag("particleFlowClusterECAL"),
    HcalPFClusters = cms.InputTag("particleFlowClusterHCAL"),
    ptMin = cms.double(1.),                     # Minimum pt                                                                         
    pMin = cms.double(1.),                      # Minimum p                                                                          
    nPixMin = cms.int32(2),                     # Nb of pixel hits                                                                   
    nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits                                                                   
    nEtaMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges                                                                
    hcalMin = cms.double(0.5),                   # Minimum hcal energy                                                               
    ecalMax = cms.double(1E9),                  # Maximum ecal energy                                                                
    verbose = cms.untracked.bool(True),         # not used.                                                                          
    #rootOutputFile = cms.string("PGun__2_200GeV__81X_upgrade2017_realistic_v22.root"),# the root tree                               
    #rootOutputFile = cms.string("PGun_8_1_0_pre16_test1.root"),# the root tree                                                       
    rootOutputFile = cms.string("step3_pi1-100_2017_realistic.root"),# the root tree                                                       
#    IsMinBias = cms.untracked.bool(False)                                                                                           
)




process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cfi")
#process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")                                                                        

process.particleFlowSimParticle.ParticleFilter = cms.PSet(
        # Allow *ALL* protons with energy > protonEMin                                                                               
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0))

process.genReReco = cms.Sequence(#process.generator+                                                                                 
                                 #process.genParticles+                                                                              
                                 #process.genJetParticles+                                                                           
                                 #process.recoGenJets+                                                                               
                                 #process.genMETParticles+                                                                           
                                 #process.recoGenMET+                                                                                
                                 process.particleFlowSimParticle)






# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction*process.hbheprerecoM3)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.bla = cms.EndPath(process.pfChargedHadronAnalyzer)
process.blo = cms.EndPath(process.genReReco)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.endjob_step,process.RECOSIMoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.endjob_step,process.blo,process.bla)

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(16)
#process.options.numberOfStreams=cms.untracked.uint32(0)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)

##from SLHCUpgradeSimulations.Configuration.HCalCustoms import load_HcalHardcode
##process = load_HcalHardcode(process)
##process.es_hardcode.useHEUpgrade = cms.bool(True)
##process.es_hardcode.useHFUpgrade = cms.bool(True)
##process.es_hardcode.heUpgrade.darkCurrent = cms.double(0)
##process.es_hardcode.SiPMCharacteristics[2].crosstalk = cms.double(0.0)
##process.es_hardcode.SiPMCharacteristics[3].crosstalk = cms.double(0.0)
##process.es_hardcode.toGet = cms.untracked.vstring('GainWidths','SiPMParameters','SiPMCharacteristics')

# Customisation from command line

dumpFile  = open("DumpRECO_Phase1_step3_GT.py", "w")
dumpFile.write(process.dumpPython())
dumpFile.close()
