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
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:./embedded_10_1_OuN.root'
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_10_1_OuN.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_1_1_Vvr.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_2_2_thx.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_3_1_lJ8.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_4_2_vF1.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_5_2_0yn.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_6_1_aYr.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_7_1_V9b.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_8_1_g2b.root',
    '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_9_2_cIx.root'
    )
    )

########################## primary vertex ###############################
#--------------------------------------------------------------------------------
# Select candidates for "the" event vertex.
#
# Selection based on PhysicsTools/PatAlogos/plugins/PATSingleVertexSelector.
# See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string
#
# NOTE: offlinePrimaryVerticesWithBS collection is sorted
#       in order of decreasing sum of Pt of tracks fitted to each vertex
#--------------------------------------------------------------------------------

# CV: ndof >= 4 if using 'offlinePrimaryVertices',
#     ndof >= 7 if using 'offlinePrimaryVerticesWithBS' as input
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
                                                    src = cms.InputTag('offlinePrimaryVertices'),
                                                    cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0"), # CV: cut >= 4 if using 'offlinePrimaryVertices',
                                                    #         >= 7 if using 'offlinePrimaryVerticesWithBS' as input
                                                    filter = cms.bool(False)                                          
                                                    )

process.selectedPrimaryVertexPosition = cms.EDFilter("VertexSelector",
                                                     src = cms.InputTag('selectedPrimaryVertexQuality'),
                                                     cut = cms.string("abs(z) < 24 & abs(position.Rho) < 2."),
                                                     filter = cms.bool(False)                                           
                                                     )

process.selectedPrimaryVertexHighestPtTrackSum = cms.EDFilter("PATSingleVertexSelector",
                                                              mode = cms.string('firstVertex'),
                                                              vertices = cms.InputTag('selectedPrimaryVertexPosition'),
                                                              filter = cms.bool(False)                                                    
                                                              )

process.selectPrimaryVertex = cms.Sequence(
    process.selectedPrimaryVertexQuality
    * process.selectedPrimaryVertexPosition
    * process.selectedPrimaryVertexHighestPtTrackSum
    )


########################## Electrons from Z ###############################
process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.load("TauAnalysis.GenSimTools.gen_decaysFromZs_cfi")


########################## analyzer ###############################

process.PFAnalyzer = cms.EDAnalyzer(
    "PFlowAnalyzer",
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    srcPFTaus = cms.InputTag("hpsPFTauProducer"),
    srcGenElectrons = cms.InputTag("genElectronsFromZs"),
    #srcGenElectrons = cms.InputTag("genElectronsFromZtautauDecays"),
    srcGenTaus = cms.InputTag("genHadronsFromZtautauDecays"),
    debug = cms.bool(False)
    )


########################## path ###############################

process.p = cms.Path(process.selectPrimaryVertex*
                     process.tauGenJets*
                     process.produceGenDecayProductsFromZs*
                     process.PFAnalyzer
                     )



########################## output ###############################
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("MVAElec.root")
    )


