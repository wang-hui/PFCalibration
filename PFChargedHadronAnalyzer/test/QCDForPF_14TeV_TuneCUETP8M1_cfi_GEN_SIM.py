# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: QCDForPF_14TeV_TuneCUETP8M1_cfi --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic50ns13TeVCollision --geometry DB:Extended --conditions 90X_upgrade2017_realistic_v15 --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('SIM',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('QCDForPF_14TeV_TuneCUETP8M1_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')

#process.generator = cms.EDFilter("Pythia8GeneratorFilter",
#    PythiaParameters = cms.PSet(
#        parameterSets = cms.vstring('pythia8CommonSettings', 
#            'pythia8CUEP8M1Settings', 
#            'processParameters'),
#        processParameters = cms.vstring('HardQCD:all = on', 
#            'PhaseSpace:pTHatMin = 15.', 
#            'PhaseSpace:pTHatMax = 3000.'),
#        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
#            'Tune:ee 7', 
#            'MultipartonInteractions:pT0Ref=2.4024', 
#            'MultipartonInteractions:ecmPow=0.25208', 
#            'MultipartonInteractions:expPow=1.6'),
#        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
#            'Main:timesAllowErrors = 10000', 
#            'Check:epTolErr = 0.01', 
#            'Beams:setProductionScalesFromLHEF = off', 
#            'SLHA:keepSM = on', 
#            'SLHA:minMassSM = 1000.', 
#            'ParticleDecays:limitTau0 = on', 
#            'ParticleDecays:tau0Max = 10', 
#            'ParticleDecays:allowPhotonRadiation = on')
#    ),
#    comEnergy = cms.double(14000.0),
#    crossSection = cms.untracked.double(74310000.0),
#    filterEfficiency = cms.untracked.double(0.037),
#    maxEventsToPrint = cms.untracked.int32(0),
#    pythiaHepMCVerbosity = cms.untracked.bool(False),
#    pythiaPylistVerbosity = cms.untracked.int32(1),
#    reweightGen = cms.PSet(

#    )
#)


process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        
        MaxE = cms.double(500),
        MaxEta = cms.double(3.0),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(200),
        MinEta = cms.double(-3.0),
        MinPhi = cms.double(-3.14159265359),
        #ParticleID = cms.vint32(211)
	PartID = cms.vint32(-211)
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single pi E 50 HCAL')
)


process.RandomNumberGeneratorService.generator.initialSeed = 15279842


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
