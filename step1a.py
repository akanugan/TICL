# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: CE_E_Front_300um_cfi -s GEN,SIM -n 10 --conditions auto:phase2_realistic_T21 --beamspot HGCALCloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry Extended2026D76 --era Phase2C11I13M9 --relval 9000,100 --no_exec --fileout file:step1.root --customise_commands from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper;randHelper=RandomNumberServiceHelper(process.RandomNumberGeneratorService);randHelper.populate()\nprocess.source.firstLuminosityBlock=cms.untracked.uint32(22)\nprocess.source.firstEvent=cms.untracked.uint32(11)
import FWCore.ParameterSet.Config as cms

import os
lumi = os.getenv('lumi')
firstEvent = os.getenv('firstEvent')

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9

from FWCore.ParameterSet.VarParsing import VarParsing

###############################################################################
### Parse input parameters ####################################################
###############################################################################
options = VarParsing('analysis')
#options = VarParsing.VarParsing('standard')
#options.register('seed',
#                 0,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Random Seed')
options.register('seed', default=None, mytype = VarParsing.varType.int)
options.register('Nevents', default=None, mytype = VarParsing.varType.int)

#options.register('filename',
#                 'OutFile.root',
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.string,
#                 'Name of the output file')
#
#options.register('Nevents',
#                 10,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Integer, number of events to be created')
#
#options.register('Nparticles',
#                 1,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Integer, number of particles per event')
#
#options.register('delta',
#                 10.0,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.float,
#                 'Float, maximum eta')
#
#options.register('overlapping',
#                 True,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.bool,
#                 'Float, maximum eta')
#
#options.register('PID',
#                 11,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Integer, PID')
#
#options.register('Emin',
#                 10.0,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.float,
#                 'Float, minimum energy in GeV')
#
#options.register('Emax',
#                 600.0,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.float,
#                 'Float, maximum energy in GeV')
#
#options.register('Etamin',
#                 1.7,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.float,
#                 'Float, minimum eta')
#
#options.register('Etamax',
#                 2.7,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.float,
#                 'Float, maximum eta')
#
options.parseArguments()

#options.register('lumi', default=None, mytype = VarParsing.varType.int)
#options.register('firstEvent', default=None, mytype = VarParsing.varType.int)

process = cms.Process('SIM',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D76_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHGCALCloseBy_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(10),
    input = cms.untracked.int32(options.Nevents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('CE_E_Front_300um_cfi nevts:10'),
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
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

'''
process.generator = cms.EDProducer("CloseByParticleGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        ControlledByEta = cms.bool(False),
        Delta = cms.double(10),
        EnMax = cms.double(200.0),
        EnMin = cms.double(25.0),
        MaxEnSpread = cms.bool(False),
        MaxEta = cms.double(2.7),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(1.7),
        MinPhi = cms.double(-3.14159265359),
        NParticles = cms.int32(1),
        Overlapping = cms.bool(False),
        PartID = cms.vint32(22),
        Pointing = cms.bool(True),
        RMax = cms.double(135.01),
        RMin = cms.double(134.99),
        RandomShoot = cms.bool(False),
        ZMax = cms.double(321.01),
        ZMin = cms.double(320.99)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('random particles in phi and r windows')
)
'''

process.generator = cms.EDProducer("CloseByParticleGunProducer",
    AddAntiParticle = cms.bool(True),
    PGunParameters = cms.PSet(
        ControlledByEta = cms.bool(True),
        Delta = cms.double(45),
        EnMax = cms.double(600),
        EnMin = cms.double(10.0),
        MaxEnSpread = cms.bool(False),
        MaxEta = cms.double(2.03),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(2.01),
        MinPhi = cms.double(-3.14159265359),
        NParticles = cms.int32(40),
        Overlapping = cms.bool(True),
        Pointing = cms.bool(True),
        PartID = cms.vint32(22,11,-11,211,-211,111,130,321,-321),
        RMax = cms.double(135.01),
        RMin = cms.double(134.99),
        RandomShoot = cms.bool(False),
        ZMax = cms.double(321.01),
        ZMin = cms.double(320.99)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('random particles in phi and r windows')
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)


process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(options.seed)
# Customisation from command line
#from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
#randHelper=RandomNumberServiceHelper(process.RandomNumberGeneratorService)
#randHelper.populate()
#process.source.firstLuminosityBlock=cms.untracked.uint32(int(lumi))
#process.source.firstEvent=cms.untracked.uint32(int(firstEvent))

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
