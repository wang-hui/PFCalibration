# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/SinglePi0E10_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO --conditions 81X_upgrade2017_realistic_v26 --era Run2_2017 --geometry Configuration.Geometry.GeometryExtended2017_cff,Configuration.Geometry.GeometryExtended2017Reco_cff --beamspot Realistic50ns13TeVCollision --eventcontent RECOSIM --customise SLHCUpgradeSimulations/Configuration/HCalCustoms.customise_Hcal2017Full -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2017_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/Generator/python/SinglePi0E10_cfi.py nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('test_phase1.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '81X_upgrade2017_realistic_v26', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxE = cms.double(200.0),
        MaxEta = cms.double(3.0),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(2.0),
        MinEta = cms.double(-3.0),
        MinPhi = cms.double(-3.14159265359),
        PartID = cms.vint32(-211)  #PDGID
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single pi0 E 10')
)

# initiate a random number seed
#process.RandomNumberGeneratorService.generator.initialSeed = xxxx
process.RandomNumberGeneratorService.generator.initialSeed = 15279842

# write variables required for calibration in a tree
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
    rootOutputFile = cms.string("PGun_TEST__2_200GeV__81X_upgrade2017_realistic_v26.root"),# the root tree
#    IsMinBias = cms.untracked.bool(False)
)

# Gen Info re-processing
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

process.genReReco = cms.Sequence(                                 
	process.particleFlowSimParticle
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.pfTree = cms.EndPath(process.pfChargedHadronAnalyzer)
process.pfSimReReco = cms.EndPath(process.genReReco)

# Schedule definition
#process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step,process.pfSimReReco,process.pfTree)
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.pfSimReReco,process.pfTree)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.HCalCustoms
from SLHCUpgradeSimulations.Configuration.HCalCustoms import customise_Hcal2017Full 

#call to customisation function customise_Hcal2017Full imported from SLHCUpgradeSimulations.Configuration.HCalCustoms
process = customise_Hcal2017Full(process)

# End of customisation functions

# Customisation from command line
