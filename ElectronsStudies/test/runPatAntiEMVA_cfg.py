from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

postfix           =  "PFlow"
runOnMC           =  True

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V14B::All')
##     process.GlobalTag.globaltag = cms.string('START42_V13::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.fileNames = cms.untracked.vstring(
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_84_2_gFP.root',
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_85_2_Qjs.root',
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_86_2_HBj.root',
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_87_2_FzR.root',
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_88_1_zDa.root'
)



################### vertex sequence ####################
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

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

################### pat specific ####################

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

from PhysicsTools.PatAlgos.tools.coreTools import *
if not runOnMC:
    removeMCMatching(process,["All"])


removeAllPATObjectsBut(process, ['Taus'])

#removeSpecificPATObjects(process, ['Photons'],
#                         outputInProcess=False)
removeCleaning(process,
               outputInProcess=False)

## #################### tau sequence #######################

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process, 
                 pfTauLabelOld = 'shrinkingConePFTauProducer',
                 pfTauLabelNew = 'hpsPFTauProducer'
                 )
process.load("RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectronMVA2_cfi")

##Define decay mode prediscriminant
requireDecayMode = cms.PSet(
    BooleanOperator = cms.string("and"),
    decayMode = cms.PSet(
        Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
        cut = cms.double(0.5)
        )    
    )
requireDecayModeAndIso = cms.PSet(
    BooleanOperator = cms.string("and"),
    decayMode = cms.PSet(
        Producer = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
        cut = cms.double(0.5)
        ),    
    Isolation = cms.PSet(
        Producer = cms.InputTag('hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr'),
        cut = cms.double(0.5)
        )
    )

process.pfRecoTauDiscriminationAgainstElectronMVA2.Prediscriminants = requireDecayMode.clone()
## process.pfRecoTauDiscriminationAgainstElectronMVA2.Prediscriminants = requireDecayModeAndIso.clone()
process.pfRecoTauDiscriminationAgainstElectronMVA2.PFTauProducer = cms.InputTag("hpsPFTauProducer")

process.makePatTaus.replace(process.patTaus,
                            process.pfRecoTauDiscriminationAgainstElectronMVA2+process.patTaus)
process.patTaus.tauIDSources= cms.PSet(
    process.patTaus.tauIDSources,
    cms.PSet(againstElectronMVA2 = cms.InputTag("pfRecoTauDiscriminationAgainstElectronMVA2") )
    )

process.selectedPatTaus.cut = "pt>15 && abs(eta)<2.3 && tauID('decayModeFinding')>0.5 && tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5"
########################## PreSelection of reconstructed Taus ###############################
process.load("RecoTauTag.TauTagTools.PFTauSelector_cfi")
process.selectedTaus = process.pfTauSelector.clone()
process.selectedTaus.src = cms.InputTag("hpsPFTauProducer")
process.selectedTaus.cut = cms.string("pt>15 && abs(eta)<2.3")
process.selectedTaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
              selectionCut=cms.double(0.5) ),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
              selectionCut=cms.double(0.5))## ,
##     cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
##               selectionCut=cms.double(0.5) )
    )
########################## Electrons and hadrons from Z ###############################
process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.load("TauAnalysis.GenSimTools.gen_decaysFromZs_cfi")

process.genElectronsPtGt10 =  cms.EDFilter("GenParticleSelector",
                                         src = cms.InputTag("genParticles"),
                                         cut = cms.string('pt > 10. & abs(pdgId) == 11'),
                                         stableOnly = cms.bool(False),
                                         filter = cms.bool(False)
                                         )

########################## Cleaned GenJets ###############################
process.genLeptonsPtGt10 =  cms.EDFilter("GenParticleSelector",
                                         src = cms.InputTag("genParticles"),
                                         cut = cms.string('pt > 10. & abs(pdgId) >= 11 & abs(pdgId) <= 16'),
                                         stableOnly = cms.bool(False),
                                         filter = cms.bool(False)
                                         )

process.genJetsAntiOverlapWithLeptonsVeto = cms.EDFilter("GenJetAntiOverlapSelector",
                                                         src = cms.InputTag('ak5GenJets'),
                                                         srcNotToBeFiltered = cms.VInputTag("genLeptonsPtGt10"),
                                                         dRmin = cms.double(0.5),
                                                         filter = cms.bool(False)
                                                         )



########################## analyzers ###############################

process.AntiEMVAAnalyzer2 = cms.EDAnalyzer(
    "AntiEMVAAnalyzer2",
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    srcPFTaus = cms.InputTag("selectedTaus"),
    srcGenElectrons = cms.InputTag("genElectronsPtGt10"),
    srcGenElectronsFromZ = cms.InputTag("genElectronsFromZs"),
    srcGenElectronsFromZTauTau = cms.InputTag("genElectronsFromZtautauDecays"),
    srcGenTaus = cms.InputTag("genHadronsFromZtautauDecays"),
    srcGenJets = cms.InputTag("genJetsAntiOverlapWithLeptonsVeto"),
    srcPatTaus = cms.InputTag("selectedPatTaus"),
    debug = cms.bool(False)
    )


process.MatchingComputerAnalyzer = cms.EDAnalyzer(
    "MatchingComputerAnalyzer",
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    srcPFTaus = cms.InputTag("selectedTaus"),
    srcGenElectrons = cms.InputTag("genElectronsPtGt10"),
    srcGenElectronsFromZ = cms.InputTag("genElectronsFromZs"),
    srcGenElectronsFromZTauTau = cms.InputTag("genElectronsFromZtautauDecays"),
    srcGenTaus = cms.InputTag("genHadronsFromZtautauDecays"),
    srcGenJets = cms.InputTag("genJetsAntiOverlapWithLeptonsVeto"),
    srcPatTaus = cms.InputTag("selectedPatTaus"),
    debug = cms.bool(False)
    )

########################## path ###############################

process.run = cms.Sequence(
    process.PFTau*
    process.selectPrimaryVertex*
    process.patDefaultSequence*
    process.tauGenJets*
    process.produceGenDecayProductsFromZs*
    process.selectedTaus*
    process.genElectronsPtGt10*
    process.genLeptonsPtGt10*process.genJetsAntiOverlapWithLeptonsVeto*
    process.AntiEMVAAnalyzer2
    
    )


process.p = cms.Path(process.run)

########################## output ###############################

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   *patEventContentNoCleaning ) 
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("AntiEMVA.root")
    )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )

process.out.fileName = cms.untracked.string('patTuplesAntiEMVA.root')

process.outpath = cms.EndPath()

## processDumpFile = open('python.dump', 'w')
## print >> processDumpFile, process.dumpPython()
