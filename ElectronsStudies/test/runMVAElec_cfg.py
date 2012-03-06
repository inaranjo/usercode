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
    #process.GlobalTag.globaltag = cms.string('START44_V5::All')
    process.GlobalTag.globaltag = cms.string('START42_V14B::All')

#    process.GlobalTag.globaltag = cms.string('START44_V10::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/user/inaranjo/DYJetsToLL-50-madgraph-tauola-Fall11-RecoTest42_PU_S6_START44_V5-v1/DYJetsToLL_46_0_4cG.root'
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_9_2_Yj1.root'
    #'/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/GEN-SIM-RECO/RecoTest44_PU_S6_START44_V5-v1/0000/FCA3A405-B9E7-E011-9BBE-E0CB4E1A1149.root'
    #'file:./embedded_10_1_OuN.root'
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_10_1_OuN.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_1_1_Vvr.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_2_2_thx.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_3_1_lJ8.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_4_2_vF1.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_5_2_0yn.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_6_1_aYr.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_7_1_V9b.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_8_1_g2b.root',
##     '/store/user/inaranjo/DYJets-50-madgraph-tauola-Chamonix12/embedded_9_2_cIx.root'
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


########################## Electrons and hadrons from Z ###############################
process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.load("TauAnalysis.GenSimTools.gen_decaysFromZs_cfi")

########################## PreSelection of reconstructed Taus ###############################
process.load("RecoTauTag.TauTagTools.PFTauSelector_cfi")
process.selectedTaus = process.pfTauSelector.clone()
process.selectedTaus.src = cms.InputTag("hpsPFTauProducer")
process.selectedTaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
              selectionCut=cms.double(0.5) ),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
              selectionCut=cms.double(0.5) )
    )

########################## Cleaned GenJets ###############################
## process.deltaRJetLeptons = cms.EDProducer(
##     "DeltaRNearestElectronComputer",
##     probes = cms.InputTag("ak5GenJets"),
##     objects = cms.InputTag("genParticles")
## ##     pdgIds = cms.vint32(-16,-15,-14,-13,-12,-11,11,12,13,14,15,16),
## ##     minPt = cms.double(10.)
##     )
## process.selectedJetsNoLeptons = cms.EDProducer(
##     "JetsCleaner",
##     jets =  cms.InputTag("ak5GenJets"),
##     valueMap = cms.InputTag("deltaRJetLeptons"),
##     minDeltaR = cms.double(0.3)
##     )

##  probes = ak5GenJets // i.e. you need to clean generator jets, while Lorenzo is cleaning reconstructed jets
##  objects = genParticles // i.e. you need to clean with respect to genParticles which have Pt > 10 GeV and are leptons (checked by pdgId)
##  pdgIds = (-16,-15,-14,-13,-12,-11,11,12,13,14,15,16)
##  minPt = cms.double(10.)

########################## analyzer ###############################

process.PFAnalyzer = cms.EDAnalyzer(
    "PFlowAnalyzer",
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    srcPFTaus = cms.InputTag("selectedTaus"),
    srcGenElectrons = cms.InputTag("genElectronsFromZs"),
    #srcGenElectrons = cms.InputTag("genElectronsFromZtautauDecays"),
    srcGenTaus = cms.InputTag("genHadronsFromZtautauDecays"),
    #srcGenJets = cms.InputTag("selectedJetsNoLeptons"),
    match = cms.int32(2),# 0 no matching, 1 match to GenElectrons, 2 match to GenTaus, 3 match to GenJets
    debug = cms.bool(False)
    )


########################## path ###############################

process.p = cms.Path(process.selectPrimaryVertex*
                     process.tauGenJets*
                     process.produceGenDecayProductsFromZs*
                     process.selectedTaus*
                     #process.deltaRJetLeptons*process.selectedJetsNoLeptons*
                     process.PFAnalyzer
                     )



########################## output ###############################
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("MVAElecTau.root")
    )


