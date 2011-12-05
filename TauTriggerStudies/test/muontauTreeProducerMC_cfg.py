import FWCore.ParameterSet.Config as cms

process = cms.Process("muontauTreeProducer")

# import of standard configurations

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# Global Tag for DoubleElectron 2011A PromptReco v1
##For DATA
##process.GlobalTag.globaltag = 'GR_P_V20::All'
##process.GlobalTag.globaltag = 'GR_P_V14::All'

#process.GlobalTag.globaltag = 'GR_E_V16::All'

##For MC
#process.GlobalTag.globaltag = 'START42_V13::All'
process.GlobalTag.globaltag = 'START42_V14B::All'


HLT_name = 'HLT'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# ---------------------------------------------------------------------
# Input File
# ---------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                            #debugFlag = cms.untracked.bool(True),
                            #debugVebosity = cms.untracked.uint32(10),
                            fileNames = cms.untracked.vstring(
    #DATA#
    #'file:/data_CMS/cms/ochando/DATA/DoubleElectron_2011A/665536EA-5E51-E011-A074-0030487C90EE.root'
    #'/store/data/Run2011A/SingleElectron/RECO/Apr22ReReco-v3/0000/0263374C-216E-E011-9948-001A92971B36.root'
    #'/store/data/Run2011A/SingleMu/AOD/May10ReReco-v1/0000/FE412F53-627B-E011-9245-0018FEFA56AE.root'
    #'/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/02DE32D4-B79C-E011-B2E4-001A4BA9B97A.root'
    '/store/mc/Summer11/DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0AB3DA9B-B3B8-E011-A81E-E0CB4E5536A8.root'
    #'/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/GEN-SIM-RECO/RecoTest44_PU_S6_START44_V5-v1/0000/FCA3A405-B9E7-E011-9BBE-E0CB4E1A1149.root'
    #'/store/data/Run2011A/SingleElectron/RECO/May10ReReco-v1/0004/404831AE-6C7B-E011-9DE2-00266CF25320.root'
    ),                         
                            )

# ---------------------------------------------------------------------
# Ouptut File
# ---------------------------------------------------------------------
process.TFileService = cms.Service ("TFileService", 
                                    fileName = cms.string ("Tree_7TeV_DATA.root")
                                    )

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

#load the EDfilter to select just skim data
#process.load("EGamma.LLRSkim.skimAllPathsFilter_cfi")
#from EGamma.LLRSkim.skimAllPathsFilter_cfi import *

#create collection of electrons with stdReco from the loose Reco
#process.load("EGamma.ECGelec.stdPreselectionSelector_cfi")
#from EGamma.ECGelec.stdPreselectionSelector_cfi import *
#process.OffRecoSelection = stdardPreselectionSelector.clone()


process.runSelection = cms.EDFilter("RunSelect",
    requireNoTimeScan = cms.untracked.bool(True) ,
    requireCollidingBX = cms.untracked.bool(False),    requireNoLumiScan = cms.untracked.bool(False),
    debug = cms.untracked.bool(False)
    )


process.recHitFilter = cms.EDFilter("RecHitFilter")

# ---------------------------------------------------------------------
# Skim ALL Path Filter
# ---------------------------------------------------------------------
#load the EDfilter to select just skim data
#process.load("EGamma.LLRSkim.skimAllPathsFilter_cfi")
#from EGamma.LLRSkim.skimAllPathsFilter_cfi import *
#process.skimAllPathsFilter = skimAllPathsFilter.clone()
#process.skimAllPathsFilter.mode = "EP"
#process.skimAllPathsFilter.mode = "ML"
#process.skimAllPathsFilter.mode = "TP"
#process.skimAllPathsFilter.eleID= "VBTF80"

# ---------------------------------------------------------------------
# JETS
# ---------------------------------------------------------------------
# JPT
#process.load('RecoJets.Configuration.RecoJPTJets_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
#JEC Corrections... to come !
# for 360: create colection of L2L3 corrected JPT jets: ak5JPTJetsL2L3  
# one need set of tags will be provided be JES
# process.p1 = cms.Path(process.ak5JPTJetsL2L3*process.dump)

# ---------------------------------------------------------------------
# Fast Jet Rho Correction
# ---------------------------------------------------------------------
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

# ---------------------------------------------------------------------
# Vertexing DA
# ---------------------------------------------------------------------
#process.load("RecoVertex.Configuration.RecoVertex_cff")
from RecoVertex.Configuration.RecoVertex_cff import *
process.vertexreco = cms.Sequence(offlinePrimaryVertices*offlinePrimaryVerticesWithBS)

# ---------------------------------------------------------------------
# Produce eID infos
# ---------------------------------------------------------------------
###process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentification_cfi")
###New optimization
#process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")

# ---------------------------------------------------------------------
# Taus
# ---------------------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")


# ---------------------------------------------------------------------
# Produce Ntuple Module
# ---------------------------------------------------------------------

produceNtuple = cms.EDAnalyzer("TriggerTnPMuonTau",
                               MuonTag = cms.InputTag("muons"),

                               VerticesTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
                               # Trigger Stuff
                               HLTTag          = cms.InputTag("TriggerResults","","HLT"),
                               HLTFilters      = cms.VInputTag("hltL1NonIsoHLTNonIsoSinglePhotonEt10HcalIsolFilter::HLT"),
                               TriggerEventTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
                               dcsTag = cms.untracked.InputTag("scalersRawToDigi"),                                 
                               type = cms.string("signalType"),
                               AOD = cms.untracked.bool(False),
                               )


process.produceNtuple = produceNtuple.clone()
process.produceNtuple.type = 'MC'
process.produceNtuple.AOD = cms.untracked.bool(False)
process.produceNtuple.FillSC = cms.untracked.bool(True)
process.produceNtuple.functionName = cms.string("EcalClusterEnergyUncertainty")
# Trigger Stuff
process.produceNtuple.HLTTag          = 'TriggerResults::' + HLT_name
process.produceNtuple.TriggerEventTag = 'hltTriggerSummaryAOD::' + HLT_name

#process.produceNtuple.HLTMuonPaths    = cms.vstring('HLT_Mu9')


process.produceNtuple.HLTFilters      = cms.VInputTag(



                                                      'hltSingleMuIsoL3IsoFiltered17::'+HLT_name, #filtre qui correspond au TAG
                                                       #ceux venant de la branche tau de mu+tau:
                                                      'hltPFTau10::'+HLT_name,
                                                      'hltPFTau10Track::'+HLT_name,
                                                      'hltFilterIsoMu12IsoPFTau10LooseIsolation::'+HLT_name,
                                                      #et au cas ou on stocke aussi ceux venant de la branche tau de tau+ele (les meme au seuil pres):
                                                      'hltPFTau15::'+HLT_name,
                                                      'hltPFTau15Track::'+HLT_name,
                                                      'hltPFTau15TrackLooseIso::'+HLT_name,
                                                      'hltPFTau20TrackLooseIso::'+HLT_name,
                                                      'hltSingleMuIsoL3IsoFiltered12::'+ HLT_name, #filtre correspondant a la partie muon du trigger mu+tau
                                                      #nouveaux filtres
                                                      'hltSingleMuIsoL3IsoFiltered15::' + HLT_name,
                                                      'hltOverlapFilterIsoMu15TightIsoPFTau20::' + HLT_name,
                                                      'hltTauJet5::' + HLT_name
                                                      
                                                      )
process.produceNtuple.HLTPaths     = cms.vstring(
                                                  'HLT_IsoMu17',#filtre qui permet de stocker les evts = TAG
                                                  'HLT_IsoMu12_LooseIsoPFTau10',
                                                  'HLT_IsoMu15_TightIsoPFTau20',#nouveau trigger
                                                  'HLT_IsoMu24',
                                                  'HLT_IsoMu24_eta2p1'
                                                 )

#process.produceNtuple.HLTFilters      = cms.VInputTag(# Muon Trigger
 #                                                     'hltSingleMu9L3Filtered9')
process.load('Htautau.TriggerStudies.pfLeptonIso_cff')


# ---------------------------------------------------------------------
# Sequence PATH
# ---------------------------------------------------------------------
process.p = cms.Path (
    process.vertexreco + # ---> to recompute Vertex
    process.kt6PFJets + ### ---> for Rho/FastJet Isolation Correction
    
    process.particleBasedIsolation +


    process.PFTau +

    process.produceNtuple 
    )



