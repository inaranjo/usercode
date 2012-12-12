import FWCore.ParameterSet.Config as cms


###### Tag ################################################################################
elePreID = "(userFloat('antiConv') > 0.5 && userFloat('nHits') <= 0 && userFloat('dxyWrtPV') < 0.02 && userFloat('dzWrtPV') < 0.1) && ((abs(superClusterPosition.Eta) < 1.479 && userFloat('sihih') < 0.01 && userFloat('dEta') < 0.007 && userFloat('dPhi') < 0.15 && userFloat('HoE') < 0.12) || (abs(superClusterPosition.Eta) > 1.479 && userFloat('sihih') < 0.03 && userFloat('dEta') < 0.009 && userFloat('dPhi') < 0.10 && userFloat('HoE') < 0.10))"

## simpleCutsVeto = "(userFloat('nHits')<=999"+ \
##                  " && (" + \
##                  " ( abs(eta)<1.5  && userFloat('sihih')<0.010 && userFloat('dPhi')<0.80 && "+ \
##                  "          userFloat('dEta') <0.007 && userFloat('HoE') <0.15)"   + \
##                  " || "  + \
##                  " ( abs(eta)>1.5 && abs(eta)<2.3 && userFloat('sihih')<0.030 && userFloat('dPhi')<0.70 && "+ \
##                  "          userFloat('dEta') <0.010 && userFloat('HoE') <999)"   + \
##                  "     )"+ \
##                  ")"
simpleCutsVeto = "userFloat('matchElectronCutsVeto')>0.5"

electronMVA = "((pt < 20 && abs(superClusterPosition.Eta) >= 0.0 && abs(superClusterPosition.Eta) < 1.0 && userFloat('mva') > 0.133) || (pt < 20 && abs(superClusterPosition.Eta) >= 1.0 && abs(superClusterPosition.Eta) < 1.5 && userFloat('mva') > 0.465) || (pt < 20 && abs(superClusterPosition.Eta) >= 1.5 && abs(superClusterPosition.Eta) < 2.5 && userFloat('mva') > 0.518) || (pt > 20 && abs(superClusterPosition.Eta) >= 0.0 && abs(superClusterPosition.Eta) < 1.0 && userFloat('mva') > 0.942) || (pt > 20 && abs(superClusterPosition.Eta) >= 1.0 && abs(superClusterPosition.Eta) < 1.5 && userFloat('mva') > 0.947) || (pt > 20 && abs(superClusterPosition.Eta) >= 1.5 && abs(superClusterPosition.Eta) < 2.5 && userFloat('mva') > 0.878))"
simpleCutsWP95 = "(userFloat('nHits')<=1  && ( (isEB && userFloat('sihih')<0.01 && userFloat('dPhi')<0.8 && userFloat('dEta')<0.007 && userFloat('HoE')<0.15) || (isEE && userFloat('sihih')<0.03 && userFloat('dPhi')<0.7 && userFloat('dEta')<0.01 && userFloat('HoE')<0.07) ))"

simpleCutsWP80 = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 &&  ( (pt>=20 && ( (isEB && userFloat('sihih')<0.01 && userFloat('dPhi')<0.06 && userFloat('dEta')<0.004 && userFloat('HoE')<0.04) || (isEE && userFloat('sihih')<0.03 && userFloat('dPhi')<0.03 && userFloat('dEta')<0.007 && userFloat('HoE')<0.025) )) || (pt<20 && (fbrem>0.15 || (abs(superClusterPosition.Eta)<1. && eSuperClusterOverP>0.95) ) && ( (isEB && userFloat('sihih')<0.01 && userFloat('dPhi')<0.03 && userFloat('dEta')<0.004 && userFloat('HoE')<0.025) || (isEE && userFloat('sihih')<0.03 && userFloat('dPhi')<0.02 && userFloat('dEta')<0.005 && userFloat('HoE')<0.025) ) )  ) )"

simpleCutsWP70 = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.03 && "  + \
                 "          userFloat('dEta')< 0.004 && userFloat('HoE') <0.025)"     + \
                 " ||"+ \
                 " (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 && " + \
                 "          userFloat('dEta') <0.005 && userFloat('HoE') <0.025) "   + \
                 "    )"   + \
                 ")"

#WW ID (HoE nell'endcap sarebbe 0.025 ma in WW e' 0.1)
eleWWID = "(userFloat('nHits')==0 && userInt('antiConv')>0.5 "+ \
                 " && ("   + \
                 " (pt>=20 && ("    + \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.06 && "  + \
                 "                        userFloat('dEta')< 0.004 && userFloat('HoE') <0.04)"     + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.007 && userFloat('HoE') <0.1) )) "+ \
                 "     || "+ \
                 " (pt<20 && (fbrem>0.15 || (abs(superClusterPosition.Eta)<1. && eSuperClusterOverP>0.95) ) && ( "+ \
                 "               (isEB && userFloat('sihih')<0.010 && userFloat('dPhi')<0.030 && " + \
                 "                        userFloat('dEta') <0.004 && userFloat('HoE') <0.025) "   + \
                 "               ||"+ \
                 "               (isEE && userFloat('sihih')<0.030 && userFloat('dPhi')<0.020 &&"  + \
                 "                        userFloat('dEta') <0.005 && userFloat('HoE') <0.1) ))" + \
                 "    )"   + \
                 ")"

tag = cms.EDFilter("PATElectronRefSelector",
                   src = cms.InputTag("selectedPatElectronsUserEmbedded"),
                   cut = cms.string('pt>25.0 && abs(eta)<2.4 && !isEBEEGap && abs(userFloat("dxyWrtPV"))<0.02 && abs(userFloat("dzWrtPV"))<0.2 && userInt("antiConv")>0.5 && userFloat("PFRelIsoDB04")<0.10'+" && "+simpleCutsWP80),
                   filter = cms.bool(False)
                   )
###### Probe ################################################################################
probe = cms.EDFilter("PATTauRefSelector",
                     src = cms.InputTag("selectedPatTausUserEmbedded"),
##                      src = cms.InputTag("patTaus"),
                     cut = cms.string('pfJetRef.pt > 20.0 && abs(eta) < 2.3'),
                     filter = cms.bool(False)
                     )

passingSecondEleVeto = probe.clone(
                        cut = cms.string(probe.cut.value() +" && "+ simpleCutsVeto+ '&& tauID("decayModeFinding")')
                        )
    
passingIsoLooseSecondEleVeto = probe.clone(
                        cut = cms.string(probe.cut.value() +" && "+ simpleCutsVeto+ '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5')
                        )
                        
                        
passingIsoLooseEleVetoLoose = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronLoose") > 0.5')
                        )

passingIsoLooseEleVetoMedium = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronMedium") > 0.5')
                        )

passingIsoLooseEleVetoTight = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronTight") > 0.5')
                        )

passingIsoLooseEleVetoMVA = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 && tauID("againstElectronMedium") > 0.5 && tauID("againstElectronMVA") > 0.5')
                        )

passingIsoMVALooseEleVetoLoose = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseIsolationMVA") > 0.5 && tauID("againstElectronLoose") > 0.5 ')
                        )

passingIsoMVALooseEleVetoMedium = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseIsolationMVA") > 0.5 && tauID("againstElectronMedium") > 0.5 ')
                        )

passingIsoMVALooseEleVetoTight = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseIsolationMVA") > 0.5 && tauID("againstElectronTight") > 0.5 ')
                        )

passingIsoMVALooseEleVetoMVA = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseIsolationMVA") > 0.5 && tauID("againstElectronMedium") > 0.5 && tauID("againstElectronMVA") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Loose1 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronLoose1MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Loose2 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronLoose2MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Medium1 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronMedium1MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Medium2 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronMedium2MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Tight1 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronTight1MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2Tight2 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronTight2MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2VTight1 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronVTight1MVA2") > 0.5')
                        )

passingIsoLooseEleVetoMVA2VTight2 = probe.clone(
                        cut = cms.string(probe.cut.value() + '&& tauID("decayModeFinding") > 0.5 && tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5  && tauID("againstElectronVTight2MVA2") > 0.5')
                        )

passingProbes = cms.Sequence(probe* 
			    (passingSecondEleVeto +
                             passingIsoLooseSecondEleVeto +
                             passingIsoLooseEleVetoLoose +
			     passingIsoLooseEleVetoMedium +
			     passingIsoLooseEleVetoTight +
			     passingIsoLooseEleVetoMVA+
			     passingIsoMVALooseEleVetoLoose+
			     passingIsoMVALooseEleVetoMedium+
			     passingIsoMVALooseEleVetoTight+
			     passingIsoMVALooseEleVetoMVA+
			     passingIsoLooseEleVetoMVA2Loose1+
			     passingIsoLooseEleVetoMVA2Loose2+
			     passingIsoLooseEleVetoMVA2Medium1+
			     passingIsoLooseEleVetoMVA2Medium2+
			     passingIsoLooseEleVetoMVA2Tight1+
			     passingIsoLooseEleVetoMVA2Tight2+
                             passingIsoLooseEleVetoMVA2VTight1+
			     passingIsoLooseEleVetoMVA2VTight2
			    )
			   )
###### Make Tag and probe pairs ################################################################################
tnp = cms.EDProducer("CandViewShallowCloneCombiner",
                     decay = cms.string("tag@+ probe@-"),
                     roles = cms.vstring('ele', 'tau'),
                     cut   = cms.string("30 < mass < 150"),
                     checkCharge = cms.bool(False),   
                     )

###### Require at least one Tag Probe pair ################################################################################
oneTp = cms.EDFilter("CandViewCountFilter",
                     src = cms.InputTag("tnp"),
                     minNumber = cms.uint32(1)
                     )
###### MC matching ################################################################################
tagMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                            pdgId = cms.vint32(11,-11),
                            src = cms.InputTag("tag"),
                            distMin = cms.double(0.15),
                            matched = cms.InputTag("genParticles")
                            )

probeMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
                            pdgId = cms.vint32(11,-11),
                            src = cms.InputTag("probe"),
                            distMin = cms.double(0.15),
                            matched = cms.InputTag("genParticles")
                            )
##### Tag and Probe Tree #################################################################################
electauTnP = cms.EDAnalyzer("TagProbeFitTreeProducer",
                            tagProbePairs = cms.InputTag("tnp"),
                            arbitration   = cms.string("OneProbe"),
                            variables = cms.PSet(
    Pt = cms.string("pt"),
    PfRelIso= cms.string('userFloat("PFRelIsoDB04")'),
    Eta =  cms.string("eta"),
    Phi = cms.string("phi"),
    DecayMode = cms.string("tauID('decayModeFinding')"),
    HpsLooseCombIsoDBCorr = cms.string("tauID('byLooseCombinedIsolationDeltaBetaCorr')"),
    HpsMediumIso = cms.string("tauID('byMediumCombinedIsolationDeltaBetaCorr')"),
    HpsTightIso = cms.string("tauID('byTightCombinedIsolationDeltaBetaCorr')"),
    HpsIsoMVARaw = cms.string("tauID('byIsolationMVAraw')"),
    HpsIsoMVALoose = cms.string("tauID('byLooseIsolationMVA')"),
    HpsIsoMVAMedium = cms.string("tauID('byMediumIsolationMVA')"),
    HpsIsoMVATight = cms.string("tauID('byTightIsolationMVA')"),
    AntiEleLoose = cms.string("tauID('againstElectronLoose')"),
    AntiEleMedium = cms.string("tauID('againstElectronMedium')"),
    AntiEleTight = cms.string("tauID('againstElectronTight')"),
    AntiEleMva = cms.string("tauID('againstElectronMVA')"),
    AntiEleMva2Loose1 = cms.string("tauID('againstElectronLoose1MVA2')"),
    AntiEleMva2Loose2 = cms.string("tauID('againstElectronLoose2MVA2')"),
    AntiEleMva2Medium1 = cms.string("tauID('againstElectronMedium1MVA2')"),
    AntiEleMva2Medium2 = cms.string("tauID('againstElectronMedium2MVA2')"),
    AntiEleMva2Tight1 = cms.string("tauID('againstElectronTight1MVA2')"),
    AntiEleMva2Tight2 = cms.string("tauID('againstElectronTight2MVA2')"),
    AntiEleMva2VTight1 = cms.string("tauID('againstElectronVTight1MVA2')"),
    AntiEleMva2VTight2 = cms.string("tauID('againstElectronVTight2MVA2')"),
    MatchElePassVeto= cms.string('userFloat("matchElectronCutsVeto")')

    ),
                            flags = cms.PSet(
    passingSecondEleVeto =  cms.InputTag("passingSecondEleVeto"),
    passingIsoLooseSecondEleVeto =  cms.InputTag("passingIsoLooseSecondEleVeto"),
    passingIsoLooseEleVetoLoose =  cms.InputTag("passingIsoLooseEleVetoLoose"),
    passingIsoLooseEleVetoMedium =  cms.InputTag("passingIsoLooseEleVetoMedium"),
    passingIsoLooseEleVetoTight =  cms.InputTag("passingIsoLooseEleVetoTight"),
    passingIsoLooseEleVetoMVA =  cms.InputTag("passingIsoLooseEleVetoMVA"),
    passingIsoMVALooseEleVetoLoose =  cms.InputTag("passingIsoMVALooseEleVetoLoose"),
    passingIsoMVALooseEleVetoMedium =  cms.InputTag("passingIsoMVALooseEleVetoMedium"),
    passingIsoMVALooseEleVetoTight =  cms.InputTag("passingIsoMVALooseEleVetoTight"),
    passingIsoMVALooseEleVetoMVA =  cms.InputTag("passingIsoMVALooseEleVetoMVA"),
    passingIsoLooseEleVetoMVA2Loose1 =  cms.InputTag("passingIsoLooseEleVetoMVA2Loose1"),
    passingIsoLooseEleVetoMVA2Loose2 =  cms.InputTag("passingIsoLooseEleVetoMVA2Loose2"),
    passingIsoLooseEleVetoMVA2Medium1 =  cms.InputTag("passingIsoLooseEleVetoMVA2Medium1"),
    passingIsoLooseEleVetoMVA2Medium2 =  cms.InputTag("passingIsoLooseEleVetoMVA2Medium2"),
    passingIsoLooseEleVetoMVA2Tight1 =  cms.InputTag("passingIsoLooseEleVetoMVA2Tight1"),
    passingIsoLooseEleVetoMVA2Tight2 =  cms.InputTag("passingIsoLooseEleVetoMVA2Tight2"),
    passingIsoLooseEleVetoMVA2VTight1 =  cms.InputTag("passingIsoLooseEleVetoMVA2VTight1"),
    passingIsoLooseEleVetoMVA2VTight2 =  cms.InputTag("passingIsoLooseEleVetoMVA2VTight2"),
    ),
                            tagVariables = cms.PSet(
    Pt              =  cms.string("pt"),                         
    Eta             =  cms.string("eta"),
    Phi             =  cms.string("phi"),
    PfRelIso        =  cms.string("userFloat('PFRelIsoDB04')"),
    GenDecay        =  cms.InputTag("addUserVariables","genDecay"),
    NumGenPU        =  cms.InputTag("addUserVariables","numGenPU"),
    ID95            =  cms.string("? (%s) ? 1.0:0.0" % simpleCutsWP95),
    ID80            =  cms.string("? (%s) ? 1.0:0.0" % simpleCutsWP80),
    ID70            =  cms.string("? (%s) ? 1.0:0.0" % simpleCutsWP70),
    EleVeto         =  cms.string("? (%s) ? 1.0:0.0" % simpleCutsVeto),
    PFRelIsoDB04    =  cms.string('userFloat("PFRelIsoDB04")'),
    AgC             =  cms.string('userInt("antiConv")'),
    HLTEle27Bit     =  cms.InputTag("addUserVariables","HLTEle27Bit"),
    HLTEle27Match   =  cms.InputTag("addUserVariables","HLTEle27Match"),
    Mt              =  cms.InputTag("addUserVariables","Mt"),
    NumEleWP70      =  cms.InputTag("addUserVariables","numEleWP70"),
    NumEleWP80      =  cms.InputTag("addUserVariables","numEleWP80"),
    NumEleWP95      =  cms.InputTag("addUserVariables","numEleWP95"),
    NumElePassVeto  =  cms.InputTag("addUserVariables","numElePassVeto"),
    PUMCWeight2011A =  cms.InputTag("addUserVariables","puMCWeightRun2011A"),
    PUMCWeight2011B =  cms.InputTag("addUserVariables","puMCWeightRun2011B"),
    PUMCWeight2012A =  cms.InputTag("addUserVariables","puMCWeightRun2012A"),

    ),
                            tagFlags = cms.PSet(),
                            pairVariables = cms.PSet(
    tnpCharge = cms.string('charge'),
    deltaR = cms.string('sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+  min( abs(daughter(0).phi-daughter(1).phi), 2*3.141 - abs(daughter(0).phi-daughter(1).phi)  ) *  min( abs(daughter(0).phi-daughter(1).phi), 2*3.141 - abs(daughter(0).phi-daughter(1).phi)  )  )')
    ),
                            pairFlags = cms.PSet(),
                            isMC = cms.bool( True ),
                            tagMatches = cms.InputTag("tagMcMatch") ,
                            probeMatches  = cms.InputTag("probeMcMatch"),
                            motherPdgId = cms.vint32(23),
                            makeMCUnbiasTree = cms.bool(True),
                            checkMotherInUnbiasEff = cms.bool(True),
                            allProbes = cms.InputTag("probe"),
                            addRunLumiInfo = cms.bool(True),
                            addEventVariablesInfo = cms.bool(True),
                         )

###### Embed user defined variables #######################################################################
vecFall11MC = (0.0020804,0.0109929,0.0228113,0.0429041,0.0566695,0.0721029,0.0729181,0.0711529,0.0713194,0.0605273,0.0531286,0.0434069,0.0446853,0.0439974,0.0422232,0.0381212,0.0286534,0.0336355,0.0270326,0.0256726,0.0236528,0.0217785,0.0180944,0.0171807,0.0145082,0.0112167,0.00852673,0.00395336,0.00458083,0.00435781,0.00354159,0.00167699,0.00103091,0.000206458,0.000618599,0.000206121,0.000412579,8.83568e-06,0,0.000412545,0,0,0,0,0,0,0,0,0,0)
vecSummer12MC = (2.344E-05,2.344E-05,2.344E-05,2.344E-05,4.687E-04,4.687E-04,7.032E-04,9.414E-04,1.234E-03,1.603E-03,2.464E-03,3.250E-03,5.021E-03,6.644E-03,8.502E-03,1.121E-02,1.518E-02,2.033E-02,2.608E-02,3.171E-02,3.667E-02,4.060E-02,4.338E-02,4.520E-02,4.641E-02,4.735E-02,4.816E-02,4.881E-02,4.917E-02,4.909E-02,4.842E-02,4.707E-02,4.501E-02,4.228E-02,3.896E-02,3.521E-02,3.118E-02,2.702E-02,2.287E-02,1.885E-02,1.508E-02,1.166E-02,8.673E-03,6.190E-03,4.222E-03,2.746E-03,1.698E-03,9.971E-04,5.549E-04,2.924E-04,1.457E-04,6.864E-05,3.054E-05,1.282E-05,5.081E-06,1.898E-06,6.688E-07,2.221E-07,6.947E-08,2.047E-08)
vecData2011A = (0,0.000115513,0.00262453,0.0231253,0.122659,0.234608,0.23849,0.171807,0.110434,0.0654693,0.0246198,0.00540239,0.000590312,5.10257e-05,2.9901e-06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
vecData2011B = (0,2.0743e-05,1.69768e-05,4.20265e-05,0.000413519,0.00233585,0.0219252,0.0625536,0.0892844,0.112663,0.130604,0.135136,0.130953,0.116884,0.0911061,0.0578308,0.0295598,0.0122804,0.00433753,0.00146337,0.000449148,0.000114735,2.54458e-05,9.34371e-07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
vecData2012A = (167025,671763,1.39495e+06,2.04506e+06,2.52639e+06,3.09467e+06,4.20416e+06,6.25511e+06,9.424e+06,1.36176e+07,1.85145e+07,2.36502e+07,2.8515e+07, 3.26439e+07,3.56836e+07,3.74294e+07,3.78322e+07,3.6981e+07,3.50694e+07,3.23539e+07,2.91133e+07,2.56141e+07,2.2085e+07,1.87038e+07,1.55924e+07,1.28217e+07,1.04196e+07,8.38235e+06,6.68523e+06,5.29173e+06,4.16073e+06,3.2513e+06,2.52558e+06,1.95021e+06,1.49674e+06,1.14137e+06,864513,650174,485349,359517,264197,192577,139220,99816.3,70974, 50050.9,35008.2,24289.3,16718.5,11417.5,7737.46,5203.99,3474.16,2302.49,1515.1,990.016,642.47,414.121,265.166,168.683,0)

addUserVariables = cms.EDProducer("UserDefinedVariablesElecTauTnP",
                                  srcElectrons = cms.InputTag("selectedPatElectronsUserEmbedded"),
                                  srcTaus = cms.InputTag("patTaus"),
                                  triggerResults = cms.InputTag("patTriggerEvent"),
                                  met = cms.InputTag("patMETsPFlow"),
                                  isMC = cms.bool(True),
                                  TrueDist2011A = cms.vdouble(vecData2011A),
                                  TrueDist2011B = cms.vdouble(vecData2011B),
                                  TrueDist2012A = cms.vdouble(vecData2012A),
                                  MCDist = cms.vdouble(vecFall11MC),
                                  debug = cms.bool(False) 
                                  )
###### Sequence #######################################################################

sequence = cms.Sequence(
    (tag+passingProbes)*
    tnp*
    oneTp+
    addUserVariables*
    (tagMcMatch+probeMcMatch)*
    electauTnP                          
    )
