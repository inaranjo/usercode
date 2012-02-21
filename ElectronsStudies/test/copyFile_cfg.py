import FWCore.ParameterSet.Config as cms

process = cms.Process("COPY")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer(
    "ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'rfio:/dpm/in2p3.fr/home/cms/trivcat/store/user/bianchi/TauPlusX/vbfSkim/fd3d8d2e69bbd9130b63d5f62133a3a0/vbfSkim_250_0_vWt.root'
    )
    )


#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:1:154411','1:1:154748','1:1:83715','1:1:80458','1:1:71441',
#    '1:1:116303','1:1:102962','1:1:2553','1:1:51804','1:1:98548',
#    '1:1:162991','1:1:31958','1:1:149135','1:1:153220','1:1:64861',
#    '1:1:127493','1:1:163868','1:1:164135','1:1:5264','1:1:178490',
#    '1:1:109878','1:1:65445','1:1:4077','1:1:30546','1:1:28180',
#    '1:1:155705','1:1:20355','1:1:76022','1:1:92733','1:1:144932',
#    '1:1:33606','1:1:147432','1:1:15632','1:1:72150','1:1:4837',
#    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

#process.p = cms.Path( process.printTree1 )

process.out = cms.OutputModule(
    "PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'keep *'),
    fileName = cms.untracked.string('embedded.root'),
)

process.outpath = cms.EndPath(process.out)
