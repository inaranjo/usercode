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
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:/data_CMS/cms/ivo/RootFiles/BJetsFall11/50D32CF1-6D3B-E111-8F8D-003048D462AE.root'
    '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_84_2_gFP.root'
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_85_2_Qjs.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_86_2_HBj.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_87_2_FzR.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_88_1_zDa.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_88_2_zAY.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_89_1_Tp2.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_89_2_Fkb.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_8_1_Ros.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_8_2_pTF.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_90_1_BOa.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_90_2_yTX.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_91_2_EvF.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_92_2_tiJ.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_93_1_xuE.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_93_2_a6f.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_94_2_fIn.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_95_2_5iA.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_96_1_KFi.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_96_2_xX4.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_98_2_bJh.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_9_2_Yj1.root',
##     '/store/user/akalinow/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/424_eletau_Fall11_v1/e8b4f85021cdba9640c984da9bbc3fb3/tautauSkimmAOD_99_2_iMt.root'
    #'/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/GEN-SIM-RECO/RecoTest44_PU_S6_START44_V5-v1/0000/FCA3A405-B9E7-E011-9BBE-E0CB4E1A1149.root'
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

process.genElectronsPtGt10 =  cms.EDFilter("GenParticleSelector",
                                         src = cms.InputTag("genParticles"),
                                         cut = cms.string('pt > 10. & abs(pdgId) == 11'),
                                         stableOnly = cms.bool(False),
                                         filter = cms.bool(False)
                                         )

########################## PreSelection of reconstructed Taus ###############################
process.load("RecoTauTag.TauTagTools.PFTauSelector_cfi")
process.selectedTaus = process.pfTauSelector.clone()
process.selectedTaus.src = cms.InputTag("hpsPFTauProducer")
process.selectedTaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
              selectionCut=cms.double(0.5) ),
##     cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),
##               selectionCut=cms.double(0.5) ),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseIsolation"),
              selectionCut=cms.double(0.5) ),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection"),
              selectionCut=cms.double(0.5) )
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
########################## GenPions ###############################
process.genPionsPtGt10 =  cms.EDFilter("GenParticleSelector",
                                         src = cms.InputTag("genParticles"),
                                         cut = cms.string('pt > 10. & abs(pdgId) == 211'),
                                         stableOnly = cms.bool(False),
                                         filter = cms.bool(False)
                                         )



########################## analyzer ###############################

process.PFAnalyzer = cms.EDAnalyzer(
    "AntiEMVAVariablesAnalyzer",
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcGsfElectrons = cms.InputTag("gsfElectrons"),
    srcPFTaus = cms.InputTag("selectedTaus"),
    srcGenElectrons = cms.InputTag("genElectronsPtGt10"),
##     srcGenElectrons = cms.InputTag("genElectronsFromZs"),
    #srcGenElectrons = cms.InputTag("genElectronsFromZtautauDecays"),
    srcGenTaus = cms.InputTag("genPionsPtGt10"),
    #srcGenTaus = cms.InputTag("genHadronsFromZtautauDecays"),
    srcGenJets = cms.InputTag("genJetsAntiOverlapWithLeptonsVeto"),
    match = cms.int32(3),# 0 no matching, 1 match to GenElectrons, 2 match to GenTaus, 3 match to GenJets
    debug = cms.bool(False)
    )


########################## path ###############################

process.p = cms.Path(process.selectPrimaryVertex*
                     process.tauGenJets*
                     process.produceGenDecayProductsFromZs*
                     process.selectedTaus*
                     process.genElectronsPtGt10*
                     process.genPionsPtGt10*
                     process.genLeptonsPtGt10*process.genJetsAntiOverlapWithLeptonsVeto*
                     process.PFAnalyzer
                     )



########################## output ###############################
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("AntiEMVAVariables.root")
    )


