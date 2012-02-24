import FWCore.ParameterSet.Config as cms

process = cms.Process("ELECMVA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

runOnMC     = True

if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START44_V10::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_9_2_cIx.root'
    )
    )
########################## analyzer ###############################

process.PFAnalyzer = cms.EDAnalyzer(
    "PFlowAnalyzer",
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    debug = cms.bool(False)
    )


########################## path ###############################

process.p = cms.Path(process.PFAnalyzer)



########################## output ###############################
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("MVAElec.root")
    )


