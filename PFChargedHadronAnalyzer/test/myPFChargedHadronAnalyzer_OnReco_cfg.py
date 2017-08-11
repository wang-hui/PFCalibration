import FWCore.ParameterSet.Config as cms

process = cms.Process("maketree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
    )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( 'file:step3.root'
        #'file:/afs/cern.ch/work/s/spandey/public/PF_cal/9_0_2_updated_PFEnergyCalibration/CMSSW_9_0_2/src/PFCalibration/PFChargedHadronAnalyzer/test/step3_RECO.root'
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseISpring17DR/SinglePion_PT200to500/AODSIM/NoPU_90X_upgrade2017_realistic_v20-v1/00000/04ED2384-5028-E711-887C-0CC47AC08C1A.root' #could copy it using xrdcp
        #'file:dipion_001758BA-CC2F-E711-927E-24BE05C6E561.root'
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_forJetMET_HCALScaleStudies', '')

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
    rootOutputFile = cms.string("step3_tree.root"),# the root tree                                                       
#    IsMinBias = cms.untracked.bool(False)                                                                                           
)


process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cfi")
process.particleFlowSimParticle.ParticleFilter = cms.PSet(
        # Allow *ALL* protons with energy > protonEMin                                                                               
        protonEMin = cms.double(5000.0),
        # Particles must have abs(eta) < etaMax (if close enough to 0,0,0)                                                           
        etaMax = cms.double(5.3),
        # Charged particles with pT < pTMin (GeV/c) are not simulated                                                                
        chargedPtMin = cms.double(0.0),
        # Particles must have energy greater than EMin [GeV]                                                                         
        EMin = cms.double(0.0))
process.genReReco = cms.Sequence(process.particleFlowSimParticle)

#process.p = cms.Path(process.demo)
process.bla = cms.Path(process.pfChargedHadronAnalyzer)
process.blo = cms.Path(process.genReReco)
process.schedule = cms.Schedule(process.blo,process.bla)

