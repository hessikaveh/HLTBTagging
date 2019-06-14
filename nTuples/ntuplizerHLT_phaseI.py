#!/usr/bin/python
"""ntuplizerHLT 
Original Code by S. Donato - https://github.com/silviodonato/usercode/tree/NtuplerFromHLT2017_V8

Code for making nTuples with offline variables (from AOD) and HLT objects (Rerun on RAW) using the heppy framework
"""
import ROOT
import itertools
import resource
import time

from array import array
from math import sqrt, pi, log10, log, exp
# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
from utils import deltaR,SetVariable,SetVariableVector,SetVariableVectorVector,DummyClass,productWithCheck,checkTriggerIndex, resetArray
from multiprocessing import Process, Queue
from copy import deepcopy

#ROOT.gROOT.LoadMacro("/scratch/sdonato/NtupleForPaolo/CMSSW_8_0_3_patch1/src/DataFormats/L1Trigger/interface/EtSumHelper.h")

Handle.productWithCheck = productWithCheck

maxJets         = 40
bunchCrossing   = 0
pt_min          = 20

def TaggedJet(jet, trackIPTagInfos ):

    small = 1.e-5
    jet_eta = jet.eta()
    jet_phi = jet.phi()

    for tag in trackIPTagInfos.productWithCheck():
        jet_p = tag.jet()
        if jet_p.isNull():
            continue;

        dr2 = (jet_p.eta()-jet_eta)**2 + (jet_p.phi()-jet_phi)**2
        if dr2 < small:
            return tag

    print "No Match Jet",jet_eta,jet_phi
    return None

def FillVectorShallowTag(shallowTag_source,variables, jetVarNames,trkVarNames,runAOD = False, offline = False, mc = False, debug = False):
    if debug: print "In FillVectorShallowTag"

    for vMember in variables.listVectorMembers:
        getattr(variables,vMember).clear()

    
    # 
    # Find the corrisponding jet info
    #
    for i in range(variables.num[0]):
        thisPhi  = variables.phi[i]
        thisEta  = variables.eta[i]
        if debug: print "thisPhi",thisPhi,"thisEta",thisEta
        
        #
        # Clear btagging variables
        #
        for jVar in jetVarNames:
            getattr(variables,jVar)[i] = -1

        for tVar in trkVarNames:
            getattr(variables,tVar).push_back(ROOT.std.vector('float')())

        if debug: print "Shallow Tag source len", len(shallowTag_source.productWithCheck())

        for obj in shallowTag_source.productWithCheck():
            tagVars = obj.taggingVariables()
        
            #
            # Get jet Pt and Eta of the btagging info
            # 
            tagJet = obj.jet()
            tagJetPhi  = tagJet.phi()
            tagJetEta  = tagJet.eta()
            #tagJetPt  = tagVars.getList(getattr(ROOT.reco.btau,"jetPt"),False)[0]
            #tagJetEta = tagVars.getList(getattr(ROOT.reco.btau,"jetEta"),False)[0]
            
            if debug: print "\ttagJetPhi",tagJetPhi,"tagJetEta",tagJetEta
            if debug: 
                for o in tagVars:
                    print "b tagging parameters:", o.first , o.second
            if (tagJetPhi == thisPhi) and (tagJetEta == thisEta):
                for jVar in jetVarNames:
                    thisJVar = tagVars.getList(getattr(ROOT.reco.btau,jVar),False)
                    if thisJVar.size() == 1:
                        getattr(variables,jVar)[i] = thisJVar[0]

                for tVar in trkVarNames:
                    thisTVar = tagVars.getList(getattr(ROOT.reco.btau,tVar),False)
                    if thisTVar.size() > 0:
                        for thisVar in thisTVar:
                            getattr(variables,tVar)[i].push_back(thisVar)


    if debug: print "Leave FillVectorShallowTag"



def FillVector(source,variables,minPt=pt_min, runAOD = False, offline = False, mc = False, debug=False):
    if debug:
        print "In FillVector before ",variables.pt
    variables.num[0] = 0
    for obj in source.productWithCheck():

        if obj.pt()<minPt:
            continue
        if variables.num[0]<len(variables.pt):
            for (name,var) in variables.__dict__.items():
                if name == "pt" :                                         var[variables.num[0]] = obj.pt()
                elif name == "eta" :                                      var[variables.num[0]] = obj.eta()
                elif name == "phi" :                                      var[variables.num[0]] = obj.phi()
                elif name == "mass" :                                     var[variables.num[0]] = obj.mass()
                elif name == "energy" :                                   var[variables.num[0]] = obj.energy()
                elif name == "neHadEF" :                                  var[variables.num[0]] = obj.neutralHadronEnergyFraction()
                elif name == "neEmEF" :                                   var[variables.num[0]] = obj.neutralEmEnergyFraction()
                elif name == "chHadEF" :                                  var[variables.num[0]] = obj.chargedHadronEnergyFraction()
                elif name == "chEmEF" :                                   var[variables.num[0]] = obj.chargedEmEnergyFraction()
                elif name == "muEF" :                                     var[variables.num[0]] = obj.muonEnergyFraction()
                elif name == "mult" :                                     var[variables.num[0]] = obj.chargedMultiplicity()+obj.neutralMultiplicity();
                elif name == "neMult" :                                   var[variables.num[0]] = obj.neutralMultiplicity()
                elif name == "chMult" :                                   var[variables.num[0]] = obj.chargedMultiplicity()
                elif name == "passesTightID":                             var[variables.num[0]] = passJetID(obj, "tight")
                elif name == "passesTightLeptVetoID":                     var[variables.num[0]] = passJetID(obj, "tightLepVeto")
                elif name == "passesLooseID":                             var[variables.num[0]] = passJetID(obj, "loose")
                elif name == 'csv' and not runAOD and offline:            var[variables.num[0]] = obj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")
                elif name == 'deepcsv_b' and not runAOD and offline:      var[variables.num[0]] = obj.bDiscriminator("pfDeepCSVJetTags:probb")
                elif name == 'deepcsv_bb' and not runAOD and offline:     var[variables.num[0]] = obj.bDiscriminator("pfDeepCSVJetTags:probbb")
                elif name == 'deepcsv_udsg' and not runAOD and offline:   var[variables.num[0]] = obj.bDiscriminator("pfDeepCSVJetTags:probudsg")
                elif name == "partonFlavour" and not runAOD and mc:       var[variables.num[0]] = obj.partonFlavour()
                elif name == "hadronFlavour" and not runAOD and mc:       var[variables.num[0]] = obj.hadronFlavour()
                
            variables.num[0] += 1

    if debug:
        print "In FillVector after ",variables.pt

def getBJetValuesforFilling(btags_source):
    bJets = {}
    if btags_source.isValid():
         bjets  = btags_source.product()
         #import pdb; pdb.set_trace()
         for ibjet in range(len(bjets)):
             bjet = bjets.key(ibjet).get()
             bJets[ibjet] = (bjet.eta(), bjet.phi(), bjets.value(ibjet))
    return bJets

def getBJetValuesforFilling_JetColl(btags_source):
    bJets = {}
    if btags_source.isValid():
         bjets  = btags_source.product()
         for ibjet in range(len(bjets)):
             bjet = bjets.key(ibjet).get()
             bJets[ibjet] = (bjet.eta(), bjet.phi(), bjets.value(ibjet))
    return bJets
    

def FillBtagAlt(bTagTuples, jets, jet_btags, jet_btagsRank = None, JetIndexVars = None, nBtagsgeNull = None, debug = False):
    """
    In this function the btags_source product is called for every time it is needed.
    For some reason, if stored (e.g. btags = btags_source.productWithCheck()), the objects
    start leaking in memory. Especially when getting the the referenced jet, this leads 
    to segmentations violations.
    """
    if debug:
        print "Running FillBtag()"
        print bTagTuples
        #print "Source is valid:",bTagTuples.isValid()
    jetB = None
    tagpairs = int(jets.num[0])*[(-1,-20)]
    for i in range(jets.num[0]):
        if debug:
            print "processing jet {0}".format(i)
        jet_btags[i] = -20.
        dRmax = 0.1
        matchcollJet = -1
        ibjet = -1
        for ibjet in range(len(bTagTuples)):
            eta, phi, tag = bTagTuples[ibjet]
            dR = deltaR(eta, phi,jets.eta[i],jets.phi[i])
            if dR<dRmax:
                jet_btags[i] = max(-1.0, tag)
                tagpairs[i] = (i, jet_btags[i])
                dRmax = dR
                matchcollJet = ibjet
        if debug:
            if dRmax < dRmax:
                print "Matched: btag coll jet {0} to online jet {1} with dR {2}".format(ibjet, i, dRmax)

    if jet_btagsRank is not None:
        if JetIndexVars is not None and isinstance(JetIndexVars, list):
            for var in JetIndexVars:
                var[0] = -1
        if nBtagsgeNull is not None:
            nBtagsgeNull[0] = 0
        from operator import itemgetter
        sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by csv value, starting with the highest
        for ipair, pair in enumerate(sortedtags):
            jet_btagsRank[pair[0]] = ipair
            if JetIndexVars is not None and isinstance(JetIndexVars, list):
                if len(JetIndexVars) >= ipair+1:
                    JetIndexVars[ipair][0] = pair[0]
            if nBtagsgeNull is not None:
                if pair[1] >= 0:
                    nBtagsgeNull[0] += 1

def FillBtagAltWLabel(bTagTuples, bTag_label, jets, jet_btags, jet_btagsRank = None, JetIndexVars = None, nBtagsgeNull = None, debug = False):
    """
    In this function the btags_source product is called for every time it is needed.
    For some reason, if stored (e.g. btags = btags_source.productWithCheck()), the objects
    start leaking in memory. Especially when getting the the referenced jet, this leads 
    to segmentations violations.
    """
    if debug:
        print "Running FillBtag()"
        print bTag_label
        print "Source is valid:",bTagTuples
    jetB = None
    tagpairs = int(jets.num[0])*[(-1,-20)]
    for i in range(jets.num[0]):
        if debug:
            print "processing jet {0}".format(i)
        jet_btags[i] = -20.
        dRmax = 0.1
        matchcollJet = -1
        ibjet = -1
        for ibjet in range(len(bTagTuples)):
            eta, phi, tag = bTagTuples[ibjet]
            dR = deltaR(eta, phi,jets.eta[i],jets.phi[i])
            if dR<dRmax:
                jet_btags[i] = max(-1.0, tag)
                tagpairs[i] = (i, jet_btags[i])
                dRmax = dR
                matchcollJet = ibjet
        if debug:
            if dRmax < dRmax:
                print "Matched: btag coll jet {0} to online jet {1} with dR {2}".format(ibjet, i, dRmax)

    if jet_btagsRank is not None:
        if JetIndexVars is not None and isinstance(JetIndexVars, list):
            for var in JetIndexVars:
                var[0] = -1
        if nBtagsgeNull is not None:
            nBtagsgeNull[0] = 0
        from operator import itemgetter
        sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by csv value, starting with the highest
        for ipair, pair in enumerate(sortedtags):
            jet_btagsRank[pair[0]] = ipair
            if JetIndexVars is not None and isinstance(JetIndexVars, list):
                if len(JetIndexVars) >= ipair+1:
                    JetIndexVars[ipair][0] = pair[0]
            if nBtagsgeNull is not None:
                if pair[1] >= 0:
                    nBtagsgeNull[0] += 1



        
def FillBtag(btags_source, jets, jet_btags, jet_btagsRank = None, JetIndexVars = None, nBtagsgeNull = None, debug = False):
    """
    In this function the btags_source product is called for every time it is needed.
    For some reason, if stored (e.g. btags = btags_source.productWithCheck()), the objects
    start leaking in memory. Especially when getting the the referenced jet, this leads 
    to segmentations violations.
    """
    if btags_source.isValid():
        if debug:
            print "Running FillBtag()"
            print "Source is valid:",btags_source.isValid()
        jetB = None
        tagpairs = int(jets.num[0])*[(-1,-20)]
        for i in range(jets.num[0]):
            if debug:
                print "processing jet {0}".format(i)
            jet_btags[i] = -20.
            dRmax = 0.1
            matchcollJet = -1
            ibjet = -1
            nbjets  = len(btags_source.product())
            for ibjet in range(nbjets):
                #jetB = btags_source.product().key(ibjet).get()
                dR = deltaR(btags_source.product().key(ibjet).get().eta(),btags_source.product().key(ibjet).get().phi(),jets.eta[i],jets.phi[i])
                #del jetB
                if dR<dRmax:
                    jet_btags[i] = max(-1.0, btags_source.product().value(ibjet))
                    tagpairs[i] = (i, jet_btags[i])
                    dRmax = dR
                    matchcollJet = ibjet
            """
            bjets = btags_source.product()
            for ibjet in range(len(bjets)):
                jetB = bjets.key(ibjet).get()
                dR = deltaR(jetB.eta(),jetB.phi(),jets.eta[i],jets.phi[i])
                del jetB
                if dR<dRmax:
                    jet_btags[i] = max(-1.0, bjets.value(ibjet))
                    tagpairs[i] = (i, jet_btags[i])
                    dRmax = dR
                    matchcollJet = ibjet
            """
            if debug:
                if dRmax < dRmax:
                    print "Matched: btag coll jet {0} to online jet {1} with dR {2}".format(ibjet, i, dRmax)

        if jet_btagsRank is not None:
            if JetIndexVars is not None and isinstance(JetIndexVars, list):
                for var in JetIndexVars:
                    var[0] = -1
            if nBtagsgeNull is not None:
                nBtagsgeNull[0] = 0
            from operator import itemgetter
            sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by csv value, starting with the highest
            for ipair, pair in enumerate(sortedtags):
                jet_btagsRank[pair[0]] = ipair
                if JetIndexVars is not None and isinstance(JetIndexVars, list):
                    if len(JetIndexVars) >= ipair+1:
                        JetIndexVars[ipair][0] = pair[0]
                if nBtagsgeNull is not None:
                    if pair[1] >= 0:
                        nBtagsgeNull[0] += 1
                    
def makeDeepCSVSumRanking(jets, variable, sumVar1, sumVar2, jet_btagsRank = None, JetIndexVars = None, nBtagsgeNull = None):
    if JetIndexVars is not None and isinstance(JetIndexVars, list):
        for var in JetIndexVars:
            var[0] = -1
    if nBtagsgeNull is not None:
        nBtagsgeNull[0] = 0
    tagpairs = int(jets.num[0])*[(-1,-20)]    
    for i in range(jets.num[0]):
        if sumVar1[i] < 0 or  sumVar2[i] < 0:
            variable[i] = -1
        else:
            variable[i] = sumVar1[i] + sumVar2[i]
        tagpairs[i] = (i, sumVar1[i] + sumVar2[i])
    from operator import itemgetter
    sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by csv value, starting with the highest
    for ipair, pair in enumerate(sortedtags):
        jet_btagsRank[pair[0]] = ipair
        if JetIndexVars is not None and isinstance(JetIndexVars, list):
            if len(JetIndexVars) >= ipair+1:
                JetIndexVars[ipair][0] = pair[0]
        if nBtagsgeNull is not None:
            if pair[1] >= 0:
                nBtagsgeNull[0] += 1

def makeCSVRanking(jets, variable, jet_btagsRank = None, JetIndexVars = None, nBtagsgeNull = None):
    if JetIndexVars is not None and isinstance(JetIndexVars, list):
        for var in JetIndexVars:
            var[0] = -1
    if nBtagsgeNull is not None:
        nBtagsgeNull[0] = 0
    tagpairs = int(jets.num[0])*[(-1,-20)]    
    for i in range(jets.num[0]):
        tagpairs[i] = (i, variable[i])
    from operator import itemgetter
    sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by csv value, starting with the highest
    for ipair, pair in enumerate(sortedtags):
        jet_btagsRank[pair[0]] = ipair
        if JetIndexVars is not None and isinstance(JetIndexVars, list):
            if len(JetIndexVars) >= ipair+1:
                JetIndexVars[ipair][0] = pair[0]
        if nBtagsgeNull is not None:
            if pair[1] >= 0:
                nBtagsgeNull[0] += 1
                
def sortJetCollection(inputcollection, outputcollection, ordervalue, saveinputorder = None):
    """
    Function to copy an an collection and reorder them according to values given as *ordervalue*.
    """
    if ordervalue not in ["csv", "deepcsv"]:
        print "new order not supported"
        return False

    #Get collection index and sortvalue
    tagpairs = int(inputcollection.num[0])*[(-1,-20)]
    for i in range(inputcollection.num[0]):
        if ordervalue == "csv":
            val = inputcollection.csv[i]
        if ordervalue == "deepcsv":
            val = inputcollection.deepcsv[i]
        tagpairs[i] = ( i, val )

    from operator import itemgetter
    sortedtags = sorted(tagpairs,key=itemgetter(1), reverse=True) #This list is ordered by (deep)csv value, starting with the highest

    outputcollection.num[0] = inputcollection.num[0]

    nonarrayvals = ["num"]
    
    for ipair, pair in enumerate(sortedtags):
        for (inputname,inputvar) in inputcollection.__dict__.items():
            for (outputname,outputvar) in outputcollection.__dict__.items():
                if inputname == outputname and inputname not in nonarrayvals:
                    #print ipair, pair[0], inputname
                    outputvar[ipair] = inputvar[pair[0]]
                    break
        if saveinputorder is not None:
            saveinputorder[ipair] = pair[0]
    
    return True

def cleanCollection(inputcollection, outputcollection, cutVariable, cutValue, indexsaveVariable, boolCut = True, verbose = False):
    passingindices = []
    for i in range(inputcollection.num[0]):
        if boolCut:
            if inputcollection.__dict__[cutVariable][i] == cutValue:
                passingindices.append(i)
        else:
            if inputcollection.__dict__[cutVariable][i] > cutValue:
                passingindices.append(i)

    nonarrayvals = ["num"]
    outputcollection.num[0] = len(passingindices)
    for iindex, index in enumerate(passingindices):
        for (inputname,inputvar) in inputcollection.__dict__.items():
            for (outputname,outputvar) in outputcollection.__dict__.items():
                if inputname == outputname and inputname not in nonarrayvals:
                    #print inputname
                    #print type(inputvar[index])
                    #print isinstance(inputvar[index], ROOT.std.vector('float'))
                    if isinstance(inputvar[index], ROOT.std.vector('float')):
                        #print "Filling",inputname
                        #print iindex
                        #print len(outputvar)
                        #print type(outputvar)
                        #outputvar[iindex] = ROOT.std.vector('float')()
                        outputvar.push_back(ROOT.std.vector('float')())
                        for vecIdx in range(inputvar[index].size()):
                            #print inputvar[index].at(vecIdx)
                            #print iindex
                            outputvar[iindex].push_back(inputvar[index].at(vecIdx))
                        #print type(inputvar[index])
                        #continue
                    outputvar[iindex] = inputvar[index]
                    break
        indexsaveVariable[iindex] = index
    if verbose:
        print "Inputcollection"
        printJetCollection(inputcollection, cutVariable)
        print "Outputcollection"
        printJetCollection(outputcollection, cutVariable)
        
    return True


def printJetCollection(inputcollection, printVar = None):
    nonarrayvals = ["num"]
    for i in range(inputcollection.num[0]):
        print "Index:",i
        if printVar is None:
            for itemname, item in inputcollection.__dict__.items():
                if itemname not in nonarrayvals:
                    print itemname,"=",item[i]
        else:
            print printVar,"=",inputcollection.__dict__[printVar][i]
    

def passJetID(jet, requestedID):
    PFJetIDLoose = False
    PFJetIDTight = False
    PFJetIDTightLepVeto = False
    if (jet.chargedMultiplicity()+jet.neutralMultiplicity()) > 1 and jet.chargedMultiplicity() > 0 and jet.chargedHadronEnergyFraction() > 0:
        if jet.neutralHadronEnergyFraction() < 0.99 and jet.neutralEmEnergyFraction() < 0.99 and jet.chargedEmEnergyFraction() < 0.99:
            PFJetIDLoose = True
        if jet.neutralHadronEnergyFraction() < 0.90 and jet.neutralEmEnergyFraction() < 0.90 and jet.chargedEmEnergyFraction() < 0.99:
            PFJetIDTight = True
            if jet.muonEnergyFraction() < 0.8 and jet.chargedEmEnergyFraction() < 0.90:
                PFJetIDTightLepVeto =  True
    if requestedID == "tight":
        return PFJetIDTight
    elif requestedID == "tightLepVeto":
        return PFJetIDTightLepVeto
    elif requestedID == "loose":
        return PFJetIDLoose

def FillMCFlavour(inputcollection, ref, refVariable, fillVariable):
    for i in range(inputcollection.num[0]):
        if ref[i] >= 0:
            fillVariable[i] = refVariable[ref[i]]
        else:
            fillVariable[i] = -99

def LeptonOverlap(jets, muons, electrons, fillVariable, DeltaR = 0.4):
    for j in range(jets.num[0]):
        overlap = False
        jetVec = ROOT.TLorentzVector()
        jetVec.SetPtEtaPhiE(jets.pt[j], jets.eta[j], jets.phi[j], jets.energy[j])
        for i in range(muons.num[0]):
            muVec = ROOT.TLorentzVector()
            muVec.SetPtEtaPhiE(muons.pt[i], muons.eta[i], muons.phi[i], muons.energy[i])

            #print jetVec.DeltaR(muVec)
            if jetVec.DeltaR(muVec) < DeltaR:
                overlap = True
                break
        if overlap is False:
            for i in range(muons.num[0]):
                elVec = ROOT.TLorentzVector()
                elVec.SetPtEtaPhiE(electrons.pt[i], electrons.eta[i], electrons.phi[i], electrons.energy[i])

                #print jetVec.DeltaR(elVec)
                if jetVec.DeltaR(elVec) < DeltaR:
                    overlap = True
                    break
        #print "Jet", j, "-",overlap
        fillVariable[j] = int(overlap)


def FillMuonVector(source, variables, vertex, muonid = "loose", debug = False):
    if vertex is None:
        return False
    variables.num[0] = 0
    if debug:
        print "Muon valid:",source.isValid()
    for obj in source.productWithCheck():
        if debug:
            print "Muon passes loose ",obj.passed(obj.CutBasedIdLoose)
            print "Muon passes medium",obj.passed(obj.CutBasedIdMedium)
            print "Muon passes tight ",obj.passed(obj.CutBasedIdTight)
        passesID = False
        if muonid == "tight":
            passesID = obj.passed(obj.CutBasedIdTight)
        if muonid == "loose":
            passesID = obj.passed(obj.CutBasedIdLoose)
        if passesID:
            for (name, var) in variables.__dict__.items():
                #See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
                if name == "pt" :           var[variables.num[0]] = obj.pt()
                elif name == "eta" :        var[variables.num[0]] = obj.eta()
                elif name == "phi" :        var[variables.num[0]] = obj.phi()
                elif name == "mass" :       var[variables.num[0]] = obj.mass()
                elif name == "energy" :     var[variables.num[0]] = obj.energy()
                elif name == "iso" :        var[variables.num[0]] = getMuonIso(obj)
                elif name == "pid" :        
                    if obj.passed(obj.CutBasedIdTight):
                        pid = 3
                    elif obj.passed(obj.CutBasedIdMedium):
                        pid = 2
                    elif obj.passed(obj.CutBasedIdLoose):
                        pid = 1
                    else: 
                        pid = 0
                    var[variables.num[0]] = pid
            variables.num[0] += 1
        return True

def getMuonIso(muon):
    MuIsoVars = muon.pfIsolationR04()
    iso = (MuIsoVars.sumChargedHadronPt + max(0.0, MuIsoVars.sumNeutralHadronEt + MuIsoVars.sumPhotonEt - 0.5 * MuIsoVars.sumPUPt)) / muon.pt();
    return iso



    
def FillElectronVector(source, variables, electronid):
    variables.num[0] = 0
    for iobj, obj in enumerate(source.productWithCheck()):
        if bool(obj.electronID(electronid)): #returns float so explicit conversion necessary
            #print obj, obj.pt(), obj.electronID(electronid)
            for (name, var) in variables.__dict__.items():
                if name == "pt" :                var[variables.num[0]] = obj.pt()
                elif name == "eta" :             var[variables.num[0]] = obj.eta()
                elif name == "phi" :             var[variables.num[0]] = obj.phi()
                elif name == "mass" :            var[variables.num[0]] = obj.mass()
                elif name == "superClusterEta" : var[variables.num[0]] = obj.superCluster().eta()
                elif name == "iso" :             var[variables.num[0]] = getElecIso(obj)
                elif name == "pid":
                    if obj.electronID("cutBasedElectronID_Fall17_94X_V1_tight"):
                        pid = 3
                    elif obj.electronID("cutBasedElectronID_Fall17_94X_V1_medium"):
                        pid = 2
                    elif obj.electronID("cutBasedElectronID_Fall17_94X_V1_loose"):
                        pid = 1
                    else: 
                        pid = 0
                    var[variables.num[0]] = pid

            variables.num[0] += 1


def getElecIso(elec):
    ElecIsoVars = elec.pfIsolationVariables();
    iso = (ElecIsoVars.sumChargedHadronPt + max(0.0, ElecIsoVars.sumNeutralHadronEt + ElecIsoVars.sumPhotonEt - 0.5 * ElecIsoVars.sumPUPt)) / elec.pt();
    return iso

        
def Matching(phi, eta, jets, debug = False):
    index = -1
    if debug:
        print "++++++++++++++++++++++++++++++++++++++"
    dRmax = 0.3
    for i in range(jets.num[0]):
        dR = deltaR(eta,phi,jets.eta[i],jets.phi[i])
        if debug:
            print i, dR
        if dR<dRmax:
            index = i
            dRmax = dR
    if debug:
        print "++++++++++++++++++++++++++++++++++++++"
    return index, dRmax

def getVertex(vertex_source):
    vertices = vertex_source.productWithCheck()
    if vertices.size()>0:
        return vertices.at(0).z()
    else:
        return -1000


def getVertices(vertex_source):
    vertices = vertex_source.productWithCheck()
    vZ0 = -1000
    nVtx = 0
    size = vertices.size()
    if size>0:
        for iVtx in range(size):
            if vertices.at(iVtx).isFake() == False and vertices.at(iVtx).ndof() > 4 and abs(vertices.at(iVtx).z()) < 24 and abs(vertices.at(iVtx).position().Rho()) < 2:
                if iVtx == 0:
                    vZ0 =  vertices.at(0).z()
                nVtx += 1
    return vZ0, nVtx

    
def WithFallback(product,method="pt"):
    if product.size()>0:
        return getattr(product[0],method)()
    else:
        return -10

def BookVector(tree,name="vector",listMembers=[],listVectorMembers=[]):
    obj = DummyClass()
    obj.num   = SetVariable(tree,name+'_num' ,'I')
    for member in listMembers:
        if "match" in member or "rank" in member or "mcFlavour" in member:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'I',name+'_num',maxJets))
        else:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'F',name+'_num',maxJets))

    obj.listVectorMembers = list(listVectorMembers)
    for vmember in listVectorMembers:
        setattr(obj,vmember,SetVariableVectorVector(tree,name+'_'+vmember))

    return obj


def getPUweight(run, truePU):
    if run not in ["RunC", "RunD", "RunE", "RunF", "RunC-F"]:
        return -1
    if run == "RunC":
        pubins = [28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
        puWeights = {28 : 1.69685876616, 29 : 1.69620554189, 30 : 1.70523385947, 31 : 1.63640795444, 32 : 1.54680792311, 33 : 1.41795946622, 34 : 1.29694884918, 35 : 1.1371673236, 36 : 0.987407342594, 37 : 0.837806841505, 38 : 0.688373995569, 39 : 0.554474841998, 40 : 0.436304816597, 41 : 0.340605346931, 42 : 0.253667659597, 43 : 0.186396020682, 44 : 0.136674652326, 45 : 0.0977525133068, 46 : 0.0680113606757, 47 : 0.0462413188393, 48 : 0.0317292824167, 49 : 0.0212675734554, 50 : 0.0141717019454, 51 : 0.00900858277463, 52 : 0.00580000207683, 53 : 0.00371004717814, 54 : 0.00229313155889, 55 : 0.00140949147272, 56 : 0.00085506576242, 57 : 0.000518662919155, 58 : 0.000311190811353, 59 : 0.00017976467166, 60 : 0.000104526213077, 61 : 6.02983285722e-05, 62 : 3.49839249797e-05}
    if run == "RunD":
        pubins = [28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
        puWeights = {28 : 2.31715121903, 29 : 2.24177384225, 30 : 2.13730741994, 31 : 1.9063698654, 32 : 1.64889913161, 33 : 1.36965123988, 34 : 1.12855498097, 35 : 0.888394784789, 36 : 0.692151282263, 37 : 0.528525685644, 38 : 0.393081277971, 39 : 0.288404482963, 40 : 0.207513300002, 41 : 0.148069939956, 42 : 0.100276686747, 43 : 0.0663950149926, 44 : 0.0433624696715, 45 : 0.0272786841925, 46 : 0.0164886759071, 47 : 0.00962849744585, 48 : 0.00561461204946, 49 : 0.0031660632755, 50 : 0.00175670796918, 51 : 0.000919612858686, 52 : 0.000481644170624, 53 : 0.000247229157998, 54 : 0.000120782284412, 55 : 5.77116414773e-05, 56 : 2.67246716448e-05, 57 : 1.21300308575e-05, 58 : 5.32878267835e-06, 59 : 2.20110441649e-06, 60 : 8.91889057314e-07, 61 : 3.48689602292e-07, 62 : 1.33062097061e-07}
    if run == "RunE":
        pubins = [28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
        puWeights = {28 : 0.996702279376, 29 : 1.02481528393, 30 : 1.08181404534, 31 : 1.11108314237, 32 : 1.14104631706, 33 : 1.14870579885, 34 : 1.16428162955, 35 : 1.14347457478, 36 : 1.128651849, 37 : 1.10914017834, 38 : 1.07806260207, 39 : 1.05015582166, 40 : 1.0213595263, 41 : 1.00637173639, 42 : 0.965094913442, 43 : 0.93054813372, 44 : 0.911157350831, 45 : 0.883919472362, 46 : 0.845033279153, 47 : 0.797091024828, 48 : 0.763257939551, 49 : 0.71540586934, 50 : 0.665655730567, 51 : 0.588341741003, 52 : 0.523327421909, 53 : 0.458988754079, 54 : 0.385837865903, 55 : 0.31994712603, 56 : 0.25984405506, 57 : 0.209515716181, 58 : 0.166026547925, 59 : 0.125922149243, 60 : 0.0955918624042, 61 : 0.071589870827, 62 : 0.0536058662532}
    if run == "RunF":
        pubins = [28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
        puWeights = {28 : 0.675357118717, 29 : 0.667386924221, 30 : 0.681069677516, 31 : 0.683691224579, 32 : 0.696125645994, 33 : 0.704787524526, 34 : 0.726815121595, 35 : 0.733422621595, 36 : 0.750407140941, 37 : 0.769102067347, 38 : 0.780229751969, 39 : 0.790320177585, 40 : 0.796661376052, 41 : 0.81597277936, 42 : 0.82349281328, 43 : 0.852923160221, 44 : 0.918617050258, 45 : 1.00064246458, 46 : 1.08821451505, 47 : 1.1721999352, 48 : 1.27608387402, 49 : 1.34520051682, 50 : 1.38684568839, 51 : 1.33473343967, 52 : 1.26850644352, 53 : 1.16608188205, 54 : 1.00825624492, 55 : 0.844582503817, 56 : 0.681422924113, 57 : 0.537840115808, 58 : 0.412145934425, 59 : 0.299583955995, 60 : 0.216879172082, 61 : 0.154866755476, 62 : 0.111179997118}
    if run == "RunC-F":
        pubins = [28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]
        puWeights = {28 : 1.22379932878, 29 : 1.21888924077, 30 : 1.22824751378, 31 : 1.19082174193, 32 : 1.14842228706, 33 : 1.08595369919, 34 : 1.03692035594, 35 : 0.962868539402, 36 : 0.901949390327, 37 : 0.844362493449, 38 : 0.784497599355, 39 : 0.73284997558, 40 : 0.686648610038, 41 : 0.65717975628, 42 : 0.620576152181, 43 : 0.600536813161, 44 : 0.603452457665, 45 : 0.614053868113, 46 : 0.626844151148, 47 : 0.638670729519, 48 : 0.663539792146, 49 : 0.673550553511, 50 : 0.67414024901, 51 : 0.634378997033, 52 : 0.593145380205, 53 : 0.539264687045, 54 : 0.463258995582, 55 : 0.387098711957, 56 : 0.312666358866, 57 : 0.24784875215, 58 : 0.191273132416, 59 : 0.140333217972, 60 : 0.102694806799, 61 : 0.0741611804725, 62 : 0.0537921537147}

    weight = -1
    #print "TruePU:",truePU
    for ibin in range(len(pubins)):
        if truePU >= pubins[ibin] and truePU < pubins[ibin+1]:
            #print "PUBin:",ibin
            #print "Central value:",pubins[ibin]
            weight = puWeights[pubins[ibin]]
    if truePU > pubins[-1]:
        weight = puWeights[pubins[-1]]

    return weight
        
##########################################################################

def launchNtupleFromHLT(fileOutput,filesInput, secondaryFiles, maxEvents,preProcessing=True, firstEvent=0, MC = False, LS = None, local = False):
    import os
    bunchCrossing   = 12
    t0 = time.time()
    print "filesInput: ",filesInput
    print "fileOutput: ",fileOutput
    print "secondaryFiles: ",secondaryFiles
    print "maxEvents: ",maxEvents
    print "LumiSections: ",LS
    print "preProcessing: ",preProcessing
    print "firstEvent: ",firstEvent
    
    doTriggerCut = True
    if doTriggerCut:
        print "+-----------------------------------------------------------------------------------------+"
        print "| TriggerCut is active. All events passing none of the triggers in the menu are discarded!|"
        print "| Note: If --setup is used only the path in the actual menu are considered for this.      |" 
        print "|       Not the ones in the setup.                                                        |"
        print "+-----------------------------------------------------------------------------------------+"
        print""
    runAOD = False
    if runAOD and False: # Don't do it!
        print "                             +----------------------------+"
        print "                             | IMPORTANT: Will run on AOD |"
        print "                             +----------------------------+"
        print""
    else:
        print "                           +--------------------------------+"
        print "                           | IMPORTANT: Will run on miniAOD |"
        print "                           +--------------------------------+"
        print""


        
    dataflags = ["MuonEG"] #NOTE: Add more flags if different data datasets are considered
        
    #isMC = bool(MC)

    if len(filesInput)>0 and (len(filter(lambda x: x in filesInput[0], dataflags)) >= 1):
        print "filesinput[0] has at least on of {0}".format(dataflags)
        isMC = False
        Signal = False
    else:
        isMC = True
        Signal = False

    print "isMC = {0}".format(isMC)

    ## Pre-processing
    if preProcessing:
        from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
        from PhysicsTools.HeppyCore.framework.config import MCComponent
        if not isMC:
            cmsRun_config = "hlt_dump_phase1.py"
        else:
            cmsRun_config = "hlt_dump_mc_phase1.py"
        if LS is not None:
            import imp

            dir_ = os.getcwd()
            cmsswConfig = imp.load_source("cmsRunProcess",os.path.expandvars(cmsRun_config))
            cmsswConfig.process.source.lumisToProcess = LS
            configfile=dir_+"/mod_"+cmsRun_config
            f = open(configfile, 'w')
            f.write(cmsswConfig.process.dumpPython())
            f.close()
            print "sed -i s/inf\)/float\( \\'inf\\'\)\)/g "+configfile
            print "sed -i s/inf,/float\(\\'inf\\'\),/g "+configfile
            os.system("sed -i s/inf\)/float\(\\'inf\\'\)\)/g "+configfile)
            os.system("sed -i s/inf,/float\(\\'inf\\'\),/g "+configfile)
            cmsRun_config = "mod_"+cmsRun_config
        print "Using: {0}".format(cmsRun_config)
        preprocessor = CmsswPreprocessor(cmsRun_config)
        cfg = MCComponent("OutputHLT",filesInput, secondaryfiles=secondaryFiles)
        print "Run cmsswPreProcessing using:"
        print cfg.name
        print cfg.files
        print cfg.secondaryfiles
        print
        try:
            preprocessor.run(cfg,".",firstEvent,maxEvents)
        except:
            print "cmsswPreProcessing failed!"
            if not local:
                print "cat cmsRun_config.py"
                config = file(cmsRun_config)
                print config.read()
                print "cat cmsRun.log"
            log = file("cmsRun.log")
            print log.read()
            #preprocessor.run(cfg,".",firstEvent,maxEvents)
            raise Exception("CMSSW preprocessor failed!")

    print "Time to preprocess: {0:10f} s".format(time.time()-t0)    
    print "Filesize of {0:8f} MB".format(os.path.getsize("cmsswPreProcessing.root") * 1e-6)

    f = ROOT.TFile(fileOutput,"recreate")
    tree = ROOT.TTree("tree","tree")

    nGenHisto = ROOT.TH1F("nGen","nGen",1,1,2)
    nPassHisto = ROOT.TH1F("nPass","nPass",1,1,2)
    
    fwLiteInputs = ["cmsswPreProcessing.root"]
    if len(filesInput)==0: exit
    import os.path
    if not os.path.isfile(fwLiteInputs[0]):
        raise Exception( fwLiteInputs[0] + " does not exist.")
    events = Events (fwLiteInputs)

    ### list of input variables ###
    ### Online
    #L1
    #l1HT_source, l1HT_label                             = Handle("BXVector<l1t::EtSum>"), ("hltGtStage2Digis","EtSum")
    #l1Jets_source, l1Jets_label                         = Handle("BXVector<l1t::Jet>"), ("hltGtStage2Digis","Jet")

    #Jets
    caloJets_source, caloJets_label                     = Handle("vector<reco::CaloJet>"), ("hltAK4CaloJetsCorrectedIDPassed") #DeepNtupler: hltAK4CaloJetsCorrected as jetToken3
    calobtag_source, calobtag_label                     = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltCombinedSecondaryVertexBJetTagsCalo")
    calodeepbtag_source, calodeepbtag_label             = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsCalo:probb")
    caloShallowTag_source, caloShallowTag_label         = Handle("vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosCalo")

    #calodeepbtag_bb_source, calodeepbtag_bb_label       = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltDeepCombinedSecondaryVertexBJetTagsCalo:probbb")
    #calodeepbtag_udsg_source, calodeepbtag_udsg_label   = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltDeepCombinedSecondaryVertexBJetTagsCalo:probudsg")
    caloPUid_source, caloPUid_label                     = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltCaloJetFromPV")

    pfJets_source, pfJets_label                         = Handle("vector<reco::PFJet>"), ("hltAK4PFJetsLooseIDCorrected") #DeepNtupler: hltPFJetForBtag as jetToken2
    #pfJets_source, pfJets_label                         = Handle("vector<reco::PFJet>"), ("hltPFJetForBtag") #DeepNtupler: hltPFJetForBtag as jetToken2
#recoJetedmRefToBaseProdTofloatsAssociationVector
    pfbtag_source, pfbtag_label                         = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltCombinedSecondaryVertexBJetTagsPF")
    pfdeepbtag_source, pfdeepbtag_label                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPF:probb")
    pfShallowTag_source, pfShallowTag_label             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfos")

    pfdeepbtag_sourcev1, pfdeepbtag_labelv1                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar1:probb")
    pfShallowTag_sourcev1, pfShallowTag_labelv1             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar1")

    pfdeepbtag_sourcev2, pfdeepbtag_labelv2                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar2:probb")
    pfShallowTag_sourcev2, pfShallowTag_labelv2             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar2")

    pfdeepbtag_sourcev3, pfdeepbtag_labelv3                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar3:probb")
    pfShallowTag_sourcev3, pfShallowTag_labelv3             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar3")

    pfdeepbtag_sourcev4, pfdeepbtag_labelv4                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar4:probb")
    pfShallowTag_sourcev4, pfShallowTag_labelv4             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar4")

    pfdeepbtag_sourcev5, pfdeepbtag_labelv5                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar5:probb")
    pfShallowTag_sourcev5, pfShallowTag_labelv5             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar5")

    pfdeepbtag_sourcev6, pfdeepbtag_labelv6                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar6:probb")
    pfShallowTag_sourcev6, pfShallowTag_labelv6             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar6")

    pfdeepbtag_sourcev7, pfdeepbtag_labelv7                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar7:probb")
    pfShallowTag_sourcev7, pfShallowTag_labelv7             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar7")

    pfdeepbtag_sourcev8, pfdeepbtag_labelv8                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar8:probb")
    pfShallowTag_sourcev8, pfShallowTag_labelv8             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar8")

    pfdeepbtag_sourcev9, pfdeepbtag_labelv9                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar9:probb")
    pfShallowTag_sourcev9, pfShallowTag_labelv9             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar9")

    pfdeepbtag_sourcev10, pfdeepbtag_labelv10                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar10:probb")
    pfShallowTag_sourcev10, pfShallowTag_labelv10             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar10")

    pfdeepbtag_sourcev11, pfdeepbtag_labelv11                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar11:probb")
    pfShallowTag_sourcev11, pfShallowTag_labelv11             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar11")

    pfdeepbtag_sourcev12, pfdeepbtag_labelv12                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar12:probb")
    pfShallowTag_sourcev12, pfShallowTag_labelv12             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar12")

    pfdeepbtag_sourcev13, pfdeepbtag_labelv13                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar13:probb")
    pfShallowTag_sourcev13, pfShallowTag_labelv13             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar13")

    pfdeepbtag_sourcev14, pfdeepbtag_labelv14                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar14:probb")
    pfShallowTag_sourcev14, pfShallowTag_labelv14             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar14")

    pfdeepbtag_sourcev15, pfdeepbtag_labelv15                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar15:probb")
    pfShallowTag_sourcev15, pfShallowTag_labelv15             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar15")

    pfdeepbtag_sourcev16, pfdeepbtag_labelv16                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar16:probb")
    pfShallowTag_sourcev16, pfShallowTag_labelv16             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar16")

    pfdeepbtag_sourcev17, pfdeepbtag_labelv17                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar17:probb")
    pfShallowTag_sourcev17, pfShallowTag_labelv17             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar17")

    pfdeepbtag_sourcev18, pfdeepbtag_labelv18                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar18:probb")
    pfShallowTag_sourcev18, pfShallowTag_labelv18             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar18")

    pfdeepbtag_sourcev19, pfdeepbtag_labelv19                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar19:probb")
    pfShallowTag_sourcev19, pfShallowTag_labelv19             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar19")



    pfdeepbtag_sourcev20, pfdeepbtag_labelv20                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar20:probb")
    pfShallowTag_sourcev20, pfShallowTag_labelv20             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar20")

    pfdeepbtag_sourcev21, pfdeepbtag_labelv21                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar21:probb")
    pfShallowTag_sourcev21, pfShallowTag_labelv21             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar21")

    pfdeepbtag_sourcev22, pfdeepbtag_labelv22                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar22:probb")
    pfShallowTag_sourcev22, pfShallowTag_labelv22             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar22")

    pfdeepbtag_sourcev23, pfdeepbtag_labelv23                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar23:probb")
    pfShallowTag_sourcev23, pfShallowTag_labelv23             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar23")

    pfdeepbtag_sourcev24, pfdeepbtag_labelv24                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar24:probb")
    pfShallowTag_sourcev24, pfShallowTag_labelv24             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar24")

    pfdeepbtag_sourcev25, pfdeepbtag_labelv25                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar25:probb")
    pfShallowTag_sourcev25, pfShallowTag_labelv25             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar25")

    pfdeepbtag_sourcev26, pfdeepbtag_labelv26                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar26:probb")
    pfShallowTag_sourcev26, pfShallowTag_labelv26             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar26")

    pfdeepbtag_sourcev27, pfdeepbtag_labelv27                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar27:probb")
    pfShallowTag_sourcev27, pfShallowTag_labelv27             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar27")

    pfdeepbtag_sourcev28, pfdeepbtag_labelv28                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar28:probb")
    pfShallowTag_sourcev28, pfShallowTag_labelv28             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar28")

    pfdeepbtag_sourcev29, pfdeepbtag_labelv29                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar29:probb")
    pfShallowTag_sourcev29, pfShallowTag_labelv29             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar29")

    pfdeepbtag_sourcev30, pfdeepbtag_labelv30                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar30:probb")
    pfShallowTag_sourcev30, pfShallowTag_labelv30             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar30")

    pfdeepbtag_sourcev31, pfdeepbtag_labelv31                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar31:probb")
    pfShallowTag_sourcev31, pfShallowTag_labelv31             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar31")

    pfdeepbtag_sourcev32, pfdeepbtag_labelv32                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar32:probb")
    pfShallowTag_sourcev32, pfShallowTag_labelv32             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar32")

    pfdeepbtag_sourcev33, pfdeepbtag_labelv33                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar33:probb")
    pfShallowTag_sourcev33, pfShallowTag_labelv33             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar33")

    pfdeepbtag_sourcev34, pfdeepbtag_labelv34                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar34:probb")
    pfShallowTag_sourcev34, pfShallowTag_labelv34             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar34")

    pfdeepbtag_sourcev35, pfdeepbtag_labelv35                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar35:probb")
    pfShallowTag_sourcev35, pfShallowTag_labelv35             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar35")

    pfdeepbtag_sourcev36, pfdeepbtag_labelv36                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar36:probb")
    pfShallowTag_sourcev36, pfShallowTag_labelv36             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar36")

    pfdeepbtag_sourcev37, pfdeepbtag_labelv37                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar37:probb")
    pfShallowTag_sourcev37, pfShallowTag_labelv37             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar37")

    pfdeepbtag_sourcev38, pfdeepbtag_labelv38                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar38:probb")
    pfShallowTag_sourcev38, pfShallowTag_labelv38             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar38")

    pfdeepbtag_sourcev39, pfdeepbtag_labelv39                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar39:probb")
    pfShallowTag_sourcev39, pfShallowTag_labelv39             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar39")

    pfdeepbtag_sourcev40, pfdeepbtag_labelv40                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar40:probb")
    pfShallowTag_sourcev40, pfShallowTag_labelv40             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar40")

    pfdeepbtag_sourcev41, pfdeepbtag_labelv41                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar41:probb")
    pfShallowTag_sourcev41, pfShallowTag_labelv41             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar41")

    pfdeepbtag_sourcev42, pfdeepbtag_labelv42                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar42:probb")
    pfShallowTag_sourcev42, pfShallowTag_labelv42             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar42")

    pfdeepbtag_sourcev43, pfdeepbtag_labelv43                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar43:probb")
    pfShallowTag_sourcev43, pfShallowTag_labelv43             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar43")

    pfdeepbtag_sourcev44, pfdeepbtag_labelv44                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar44:probb")
    pfShallowTag_sourcev44, pfShallowTag_labelv44             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar44")

    pfdeepbtag_sourcev45, pfdeepbtag_labelv45                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar45:probb")
    pfShallowTag_sourcev45, pfShallowTag_labelv45             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar45")

    pfdeepbtag_sourcev46, pfdeepbtag_labelv46                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar46:probb")
    pfShallowTag_sourcev46, pfShallowTag_labelv46             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar46")

    pfdeepbtag_sourcev47, pfdeepbtag_labelv47                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar47:probb")
    pfShallowTag_sourcev47, pfShallowTag_labelv47             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar47")

    pfdeepbtag_sourcev48, pfdeepbtag_labelv48                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar48:probb")
    pfShallowTag_sourcev48, pfShallowTag_labelv48             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar48")

    pfdeepbtag_sourcev49, pfdeepbtag_labelv49                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar49:probb")
    pfShallowTag_sourcev49, pfShallowTag_labelv49             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar49")

    pfdeepbtag_sourcev50, pfdeepbtag_labelv50                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar50:probb")
    pfShallowTag_sourcev50, pfShallowTag_labelv50             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar50")

    pfdeepbtag_sourcev51, pfdeepbtag_labelv51                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar51:probb")
    pfShallowTag_sourcev51, pfShallowTag_labelv51             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar51")

    pfdeepbtag_sourcev52, pfdeepbtag_labelv52                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar52:probb")
    pfShallowTag_sourcev52, pfShallowTag_labelv52             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar52")

    pfdeepbtag_sourcev53, pfdeepbtag_labelv53                 = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>>"), ("hltDeepCombinedSecondaryVertexBJetTagsPFvar53:probb")
    pfShallowTag_sourcev53, pfShallowTag_labelv53             = Handle("std::vector<reco::ShallowTagInfo>"),("hltDeepCombinedSecondaryVertexBJetTagsInfosvar53")



    #pfdeepbtag_bb_source, pfdeepbtag_bb_label           = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltDeepCombinedSecondaryVertexBJetTagsPF:probbb")
    #pfdeepbtag_udsg_source, pfdeepbtag_udsg_label       = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltDeepCombinedSecondaryVertexBJetTagsPF:probudsg")

    #MET and HT
    #pfMet_source, pfMet_label                           = Handle("vector<reco::PFMET>"), ("hltPFMETProducer") #Not working. .product() throws expection
    #pfMht_source, pfMht_label                           = Handle("vector<reco::MET>"), ("hltPFMHTTightID") #Not working. .product() throws expection
    #caloMet_source, caloMet_label                       = Handle("vector<reco::CaloMET>"), ("hltMet") #Not working. .product() throws expection
    #caloMht_source, caloMht_label                       = Handle("vector<reco::MET>"), ("hltMht") #Not working. .product() throws expection
    #caloMhtNoPU_source, caloMhtNoPU_label               = Handle("vector<reco::MET>"), ("hltMHTNoPU") #Not working. .product() throws expection

    #Vertex
    FastPrimaryVertex_source, FastPrimaryVertex_label   = Handle("vector<reco::Vertex>"), ("hltFastPrimaryVertex")
    FPVPixelVertices_source, FPVPixelVertices_label     = Handle("vector<reco::Vertex>"), ("hltFastPVPixelVertices")
    PixelVertices_source, PixelVertices_label           = Handle("vector<reco::Vertex>"), ("hltPixelVertices")
    VerticesPF_source, VerticesPF_label                 = Handle("vector<reco::Vertex>"), ("hltVerticesPF")
    VerticesL3_source, VerticesL3_label                 = Handle("vector<reco::Vertex>"), ("hltVerticesL3")

    #The rest
    #triggerBits, triggerBitLabel                        = Handle("edm::TriggerResults"), ("TriggerResults::MYHLT")
    triggerBits, triggerBitLabel                        = Handle("edm::TriggerResults"), ("TriggerResults::HLT")



    if runAOD:
        #Leptons
        offEle_source, offEle_label                         = Handle("vector<reco::GsfElectron>"), ("gedGsfElectrons")
        offMu_source, offMu_label                           = Handle("vector<reco::Muon>"), ("muons")
        MuGlobalTracks_source, MuGlobalTracks_label         = Handle("vector<reco::Track>"), ("globalTracks")

        #Jets
        offJets_source, offJets_label                       = Handle("vector<reco::PFJet>"), ("ak4PFJetsCHS") #DeepNtuple: jetToken1
        #offJetsnoCHS_source, offJetsnoCHS_label             = Handle("vector<reco::PFJet>"), ("ak4PFJets")
        offbtag_source, offbtag_label                       = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        offdeepbtag_source, offdeepbtag_label               = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("pfDeepCSVJetTags:probb")
        offdeepbtag_bb_source, offdeepbtag_bb_label         = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("pfDeepCSVJetTags:probbb")
        offdeepbtag_udsg_source, offdeepbtag_udsg_label     = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("pfDeepCSVJetTags:probudsg")

        #MET and HT
        offMet_source, offMet_label                         = Handle("vector<reco::PFMET>"), ("pfMet")

        #Vertex
        pileUp_source, pileUp_label                         = Handle("vector<PileupSummaryInfo>"), ("addPileupInfo")
        VerticesOff_source, VerticesOff_label               = Handle("vector<reco::Vertex>"), ("offlinePrimaryVertices")

        #Gen
        #genJets_source, genJets_label                       = Handle("vector<reco::GenJet>"), ("ak4GenJetsNoNu")
        #genMet_source, genMet_label                         = Handle("vector<reco::GenMET>"), ("genMetTrue")
        genParticles_source, genParticles_label             = Handle("vector<reco::GenParticle>"), ("genParticles")
        generator_source, generator_label                   = Handle("GenEventInfoProduct"), ("generator")
    else:
        if isMC:
            offEle_source, offEle_label                         = Handle("vector<pat::Electron>"), ("slimmedElectrons::RECO") #NOTE: This should match the process name in your hlt_dump!
            idName = "cutBasedElectronID_Fall17_94X_V1_loose" #in MINIAODv2 change _ to -
        else:
            offEle_source, offEle_label                         = Handle("vector<pat::Electron>"), ("slimmedElectrons::MYHLT")
            idName = "cutBasedElectronID_Fall17_94X_V1_loose" #in MINIAODv2 change _ to -
            #cutBasedElectronID_Fall17_94X_V1_loose
            #cutBasedElectronID_Fall17_94X_V1_medium 
            #cutBasedElectronID_Fall17_94X_V1_tight

        offMu_source, offMu_label                           = Handle("vector<pat::Muon>"), ("slimmedMuons")
        MuGlobalTracks_source, MuGlobalTracks_label         = Handle("vector<reco::Track>"), ("globalTracks")

        #Jets
        offJets_source, offJets_label                       = Handle("vector<pat::Jet>"), ("slimmedJets")
        offShallowTag_source, offShallowTag_label           = Handle("vector<reco::ShallowTagInfo>"), ("TESTpfDeepCSVTagInfos")

        #MET
        offMet_source, offMet_label                         = Handle("std::vector<pat::MET>"), ("slimmedMETs")

        #Vertex
        pileUp_source, pileUp_label                         = Handle("vector<PileupSummaryInfo>"), ("slimmedAddPileupInfo")
        VerticesOff_source, VerticesOff_label               = Handle("vector<reco::Vertex>"), ("offlineSlimmedPrimaryVertices")

        #Gen
        #genJets_source, genJets_label                       = Handle("vector<reco::GenJet>"), ("slimmedGenJets")
        #genMet_source, genMet_label                         = Handle("vector<reco::GenMET>"), ("genMetTrue")
        genParticles_source, genParticles_label             = Handle("vector<reco::GenParticle>"), ("prunedGenParticles")
        generator_source, generator_label                   = Handle("GenEventInfoProduct"), ("generator")


    jetVars = ["jetEnergy", "jetPt","jetEta","jetAbsEta","jetPhi","jetNTracks","jetNTracksEtaRel",
               "jetNSecondaryVertices","jetNSelectedTracks","vertexMass","vertexNTracks","vertexFitProb",
               "vertexEnergyRatio", "vertexJetDeltaR", "trackJetPt","trackSumJetEtRatio","trackSumJetDeltaR",         
               "vertexCategory",
               "vertexLeptonCategory",
               "leptonQuality",
               "leptonQuality2",
               "chargedHadronEnergyFraction",
               "neutralHadronEnergyFraction",
               "photonEnergyFraction",
               "electronEnergyFraction",
               "muonEnergyFraction",
               "chargedHadronMultiplicity",
               "neutralHadronMultiplicity",
               "photonMultiplicity",
               "electronMultiplicity",
               "muonMultiplicity",
               "hadronMultiplicity",
               "hadronPhotonMultiplicity",
               "massVertexEnergyFraction",
               "totalMultiplicity",
               "vertexBoostOverSqrtJetPt",
               "leptonSip2d",
               "leptonSip3d",
               "leptonPtRel",
               "leptonP0Par",
               "leptonEtaRel",
               "leptonDeltaR",
               "leptonRatio",
               "leptonRatioRel",
               "electronMVA",
               "Jet_SoftMu",
               "Jet_SoftEl",
               "Jet_JBP",
               "Jet_JP",
               "flightDistance1dVal","flightDistance1dSig",
               "flightDistance2dVal","flightDistance2dSig",
               "flightDistance3dVal","flightDistance3dSig",
               "trackSip2dValAboveCharm","trackSip2dSigAboveCharm",
               "trackSip3dValAboveCharm","trackSip3dSigAboveCharm",
               ]


    trkVars = ["trackSip3dVal",
               "trackSip3dSig",
               "trackSip2dVal",
               "trackSip2dSig",
               "trackDecayLenVal",
               "trackDecayLenSig",
               "trackJetDistVal",
               "trackJetDistSig",
               "trackGhostTrackWeight",
               "trackGhostTrackDistSig",
               "trackGhostTrackDistVal",
               "trackPtRel",
               "trackMomentum",
               "trackEta",
               "trackPhi",
               "trackCharge",
               "trackPPar",
               "trackDeltaR",
               "trackEtaRel",
               "trackPtRatio",
               "trackPParRatio",
               "trackP0Par",
               "trackP0ParRatio",
               "trackChi2",
               "trackNTotalHits",
               "trackNPixelHits",
               ]
        


    ### create output variables ###
    #Leptons
    offTightElectrons   = BookVector(tree, "offTightElectrons", ['pt','eta', 'phi','mass', 'energy', "superClusterEta",'iso','pid'])
    offTightMuons       = BookVector(tree, "offTightMuons", ['pt','eta', 'phi','mass', 'energy', 'iso','pid'])
    
    #Jets:
    #l1Jets              = BookVector(tree,"l1Jets",['pt','eta','phi','energy','matchOff','matchGen'])
    caloJets            = BookVector(tree,"caloJets",['pt','eta','phi','mass', 'energy','matchOff','matchGen','puId','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJets              = BookVector(tree,"pfJets",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv1              = BookVector(tree,"pfJetsv1",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv2              = BookVector(tree,"pfJetsv2",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv3              = BookVector(tree,"pfJetsv3",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv4              = BookVector(tree,"pfJetsv4",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv5              = BookVector(tree,"pfJetsv5",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv6              = BookVector(tree,"pfJetsv6",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv7              = BookVector(tree,"pfJetsv7",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv8              = BookVector(tree,"pfJetsv8",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv9              = BookVector(tree,"pfJetsv9",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv10              = BookVector(tree,"pfJetsv10",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)


    pfJetsv11              = BookVector(tree,"pfJetsv11",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv12              = BookVector(tree,"pfJetsv12",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv13              = BookVector(tree,"pfJetsv13",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv14              = BookVector(tree,"pfJetsv14",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv15              = BookVector(tree,"pfJetsv15",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv16              = BookVector(tree,"pfJetsv16",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv17              = BookVector(tree,"pfJetsv17",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv18              = BookVector(tree,"pfJetsv18",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv19              = BookVector(tree,"pfJetsv19",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv20              = BookVector(tree,"pfJetsv20",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv21              = BookVector(tree,"pfJetsv21",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv22              = BookVector(tree,"pfJetsv22",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv23              = BookVector(tree,"pfJetsv23",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv24              = BookVector(tree,"pfJetsv24",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv25              = BookVector(tree,"pfJetsv25",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv26              = BookVector(tree,"pfJetsv26",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv27              = BookVector(tree,"pfJetsv27",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv28              = BookVector(tree,"pfJetsv28",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv29              = BookVector(tree,"pfJetsv29",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv30              = BookVector(tree,"pfJetsv30",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv31              = BookVector(tree,"pfJetsv31",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv32              = BookVector(tree,"pfJetsv32",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv33              = BookVector(tree,"pfJetsv33",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv34              = BookVector(tree,"pfJetsv34",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv35              = BookVector(tree,"pfJetsv35",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv36              = BookVector(tree,"pfJetsv36",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv37              = BookVector(tree,"pfJetsv37",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv38              = BookVector(tree,"pfJetsv38",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv39              = BookVector(tree,"pfJetsv39",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv40              = BookVector(tree,"pfJetsv40",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv41              = BookVector(tree,"pfJetsv41",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)
    pfJetsv42              = BookVector(tree,"pfJetsv42",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv43              = BookVector(tree,"pfJetsv43",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv44              = BookVector(tree,"pfJetsv44",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv45              = BookVector(tree,"pfJetsv45",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv46              = BookVector(tree,"pfJetsv46",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv47              = BookVector(tree,"pfJetsv47",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv48              = BookVector(tree,"pfJetsv48",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv49              = BookVector(tree,"pfJetsv49",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv50              = BookVector(tree,"pfJetsv50",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv51              = BookVector(tree,"pfJetsv51",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv52              = BookVector(tree,"pfJetsv52",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)

    pfJetsv53              = BookVector(tree,"pfJetsv53",['pt','eta','phi','mass', 'energy','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','csv','deepcsv','deepcsv_bb','deepcsv_udsg',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "mcFlavour"]+jetVars,trkVars)


    offJets             = BookVector(tree,"offJets",['pt','eta','phi','mass', 'energy','csv','deepcsv','deepcsv_bb','deepcsv_b','deepcsv_udsg','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "matchPF", "matchCalo", "mcFlavour", "partonFlavour", "hadronFlavour", "lepOverlap04Tight"]+jetVars,trkVars)
    #offCleanJets        = BookVector(tree,"offCleanJets",['pt','eta','phi','mass', 'energy','csv','deepcsv','deepcsv_bb','deepcsv_b','deepcsv_udsg','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankCSV", "rankDeepCSV", "matchPF", "matchCalo", "mcFlavour", "partonFlavour", "hadronFlavour", "lepOverlap04Tight", "offOrder"]+jetVars,trkVars)
    #offCleanCSVJets          = BookVector(tree,"offCleanCSVJets",['pt','eta','phi','mass', 'energy','csv','deepcsv','deepcsv_bb','deepcsv_b','deepcsv_udsg','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankpt", "matchPF", "matchCalo", "mcFlavour", "partonFlavour", "hadronFlavour", "lepOverlap04Tight"])
    #offCleanDeepCSVJets      = BookVector(tree,"offCleanDeepCSVJets",['pt','eta','phi','mass', 'energy','csv','deepcsv','deepcsv_bb','deepcsv_b','deepcsv_udsg','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult',"passesTightID","passesTightLeptVetoID", "passesLooseID", "rankpt", "matchPF", "matchCalo", "mcFlavour", "partonFlavour", "hadronFlavour", "lepOverlap04Tight"])

    """
    if isMC:
        genJets             = BookVector(tree,"genJets",['pt','eta','phi','mass','mcFlavour','mcPt'])
    """
    #MET and HT
    #l1HT                = SetVariable(tree,'l1HT')
    caloMet             = BookVector(tree,"caloMet",['pt','phi'])
    caloMht             = BookVector(tree,"caloMht",['pt','phi'])
    caloMhtNoPU         = BookVector(tree,"caloMhtNoPU",['pt','phi'])
    pfMet               = BookVector(tree,"pfMet",['pt','phi'])
    pfMht               = BookVector(tree,"pfMht",['pt','phi'])
    l1Met               = SetVariable(tree,'l1Met')
    l1Met_phi           = SetVariable(tree,'l1Met_phi')
    l1Mht               = SetVariable(tree,'l1Mht')
    l1Mht_phi           = SetVariable(tree,'l1Mht')
    offMet              = BookVector(tree,"offMet",['pt','phi'])
    if isMC:
        genMet              = BookVector(tree,"genMet",['pt','phi'])
    
    #Vertex
    FastPrimaryVertex   = SetVariable(tree,'FastPrimaryVertex')
    FPVPixelVertices    = SetVariable(tree,'FPVPixelVertices')
    PixelVertices       = SetVariable(tree,'PixelVertices')
    VerticesPF          = SetVariable(tree,'VerticesPF')
    VerticesL3          = SetVariable(tree,'VerticesL3')
    VerticesOff         = SetVariable(tree,'VerticesOff')
    nOffVertices        = SetVariable(tree,"nPV")
    trueVertex          = SetVariable(tree,'trueVertex')
    
    #General event variables
    evt                 = SetVariable(tree,'evt')
    lumi                = SetVariable(tree,'lumi')
    run                 = SetVariable(tree,'run')

    if isMC:
        pu              = SetVariable(tree,'pu')
        ptHat           = SetVariable(tree,'ptHat')
        maxPUptHat      = SetVariable(tree,'maxPUptHat')
        #PU weights
        wPURunC         = SetVariable(tree, "wPURunC")
        wPURunD         = SetVariable(tree, "wPURunD")
        wPURunE         = SetVariable(tree, "wPURunE")
        wPURunF         = SetVariable(tree, "wPURunF")
        wPURunCF        = SetVariable(tree, "wPURunCF")


    f.cd()
    crun = 0
    cls = 0
    ##get trigger names
    events.to(0)
    print events
    for event in events: break
    event.getByLabel(triggerBitLabel, triggerBits)
    names = event.object().triggerNames(triggerBits.product())
    triggerNames = names.triggerNames()
    for name in triggerNames: name = name.split("_v")[0]
    nTriggers = len(triggerNames)
    triggerVars = {}
    for trigger in triggerNames:
        if trigger.startswith("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_") or trigger.startswith("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_"):
            triggerVars[trigger]=array( 'i', [ 0 ] )
            tree.Branch( trigger, triggerVars[trigger], trigger+'/O' )

    ##event loop
    print "Starting event loop"
    for iev,event in enumerate(events):
        #raw_input("start event")
        if iev>maxEvents and maxEvents>=0: break
        #print "Starting with event",iev
        nGenHisto.Fill(1)
        #print "Event: {0}".format(iev)
        ####################################################
        ####################################################
        #Getting L1 handles
        #event.getByLabel(l1Jets_label, l1Jets_source)
        #event.getByLabel(l1HT_label, l1HT_source)

        #Getting Lepton handles
        event.getByLabel(offEle_label, offEle_source)
        event.getByLabel(offMu_label, offMu_source)
        event.getByLabel(MuGlobalTracks_label, MuGlobalTracks_source)

        #Getting Jet handles
        event.getByLabel(caloJets_label, caloJets_source)
        event.getByLabel(calobtag_label, calobtag_source)
        event.getByLabel(calodeepbtag_label, calodeepbtag_source)
        #event.getByLabel(calodeepbtag_bb_label, calodeepbtag_bb_source)
        #event.getByLabel(calodeepbtag_udsg_label, calodeepbtag_udsg_source)
        event.getByLabel(caloPUid_label, caloPUid_source)
        event.getByLabel(caloShallowTag_label, caloShallowTag_source)
        
        event.getByLabel(pfJets_label, pfJets_source)
        event.getByLabel(pfbtag_label, pfbtag_source)
        event.getByLabel(pfdeepbtag_label, pfdeepbtag_source)
        #event.getByLabel(pfdeepbtag_bb_label, pfdeepbtag_bb_source)
        #event.getByLabel(pfdeepbtag_udsg_label, pfdeepbtag_udsg_source)
        event.getByLabel(pfShallowTag_label, pfShallowTag_source)

        event.getByLabel(pfdeepbtag_labelv1, pfdeepbtag_sourcev1)
        event.getByLabel(pfShallowTag_labelv1, pfShallowTag_sourcev1)
        event.getByLabel(pfdeepbtag_labelv2, pfdeepbtag_sourcev2)
        event.getByLabel(pfShallowTag_labelv2, pfShallowTag_sourcev2)
        event.getByLabel(pfdeepbtag_labelv3, pfdeepbtag_sourcev3)
        event.getByLabel(pfShallowTag_labelv3, pfShallowTag_sourcev3)
        event.getByLabel(pfdeepbtag_labelv4, pfdeepbtag_sourcev4)
        event.getByLabel(pfShallowTag_labelv4, pfShallowTag_sourcev4)
        event.getByLabel(pfdeepbtag_labelv5, pfdeepbtag_sourcev5)
        event.getByLabel(pfShallowTag_labelv5, pfShallowTag_sourcev5)
        event.getByLabel(pfdeepbtag_labelv6, pfdeepbtag_sourcev6)
        event.getByLabel(pfShallowTag_labelv6, pfShallowTag_sourcev6)
        event.getByLabel(pfdeepbtag_labelv7, pfdeepbtag_sourcev7)
        event.getByLabel(pfShallowTag_labelv7, pfShallowTag_sourcev7)
        event.getByLabel(pfdeepbtag_labelv8, pfdeepbtag_sourcev8)
        event.getByLabel(pfShallowTag_labelv8, pfShallowTag_sourcev8)
        event.getByLabel(pfdeepbtag_labelv9, pfdeepbtag_sourcev9)
        event.getByLabel(pfShallowTag_labelv9, pfShallowTag_sourcev9)
        event.getByLabel(pfdeepbtag_labelv10, pfdeepbtag_sourcev10)
        event.getByLabel(pfShallowTag_labelv10, pfShallowTag_sourcev10)

        event.getByLabel(pfdeepbtag_labelv11, pfdeepbtag_sourcev11)
        event.getByLabel(pfShallowTag_labelv11, pfShallowTag_sourcev11)
        event.getByLabel(pfdeepbtag_labelv12, pfdeepbtag_sourcev12)
        event.getByLabel(pfShallowTag_labelv12, pfShallowTag_sourcev12)
        event.getByLabel(pfdeepbtag_labelv13, pfdeepbtag_sourcev13)
        event.getByLabel(pfShallowTag_labelv13, pfShallowTag_sourcev13)
        event.getByLabel(pfdeepbtag_labelv14, pfdeepbtag_sourcev14)
        event.getByLabel(pfShallowTag_labelv14, pfShallowTag_sourcev14)
        event.getByLabel(pfdeepbtag_labelv15, pfdeepbtag_sourcev15)
        event.getByLabel(pfShallowTag_labelv15, pfShallowTag_sourcev15)
        event.getByLabel(pfdeepbtag_labelv16, pfdeepbtag_sourcev16)
        event.getByLabel(pfShallowTag_labelv16, pfShallowTag_sourcev16)
        event.getByLabel(pfdeepbtag_labelv17, pfdeepbtag_sourcev17)
        event.getByLabel(pfShallowTag_labelv17, pfShallowTag_sourcev17)
        event.getByLabel(pfdeepbtag_labelv18, pfdeepbtag_sourcev18)
        event.getByLabel(pfShallowTag_labelv18, pfShallowTag_sourcev18)
        event.getByLabel(pfdeepbtag_labelv19, pfdeepbtag_sourcev19)
        event.getByLabel(pfShallowTag_labelv19, pfShallowTag_sourcev19)
        event.getByLabel(pfdeepbtag_labelv20, pfdeepbtag_sourcev20)
        event.getByLabel(pfShallowTag_labelv20, pfShallowTag_sourcev20)

        event.getByLabel(pfdeepbtag_labelv21, pfdeepbtag_sourcev21)
        event.getByLabel(pfShallowTag_labelv21, pfShallowTag_sourcev21)
        event.getByLabel(pfdeepbtag_labelv22, pfdeepbtag_sourcev22)
        event.getByLabel(pfShallowTag_labelv22, pfShallowTag_sourcev22)
        event.getByLabel(pfdeepbtag_labelv23, pfdeepbtag_sourcev23)
        event.getByLabel(pfShallowTag_labelv23, pfShallowTag_sourcev23)
        event.getByLabel(pfdeepbtag_labelv24, pfdeepbtag_sourcev24)
        event.getByLabel(pfShallowTag_labelv24, pfShallowTag_sourcev24)
        event.getByLabel(pfdeepbtag_labelv25, pfdeepbtag_sourcev25)
        event.getByLabel(pfShallowTag_labelv25, pfShallowTag_sourcev25)
        event.getByLabel(pfdeepbtag_labelv26, pfdeepbtag_sourcev26)
        event.getByLabel(pfShallowTag_labelv26, pfShallowTag_sourcev26)
        event.getByLabel(pfdeepbtag_labelv27, pfdeepbtag_sourcev27)
        event.getByLabel(pfShallowTag_labelv27, pfShallowTag_sourcev27)
        event.getByLabel(pfdeepbtag_labelv28, pfdeepbtag_sourcev28)
        event.getByLabel(pfShallowTag_labelv28, pfShallowTag_sourcev28)
        event.getByLabel(pfdeepbtag_labelv29, pfdeepbtag_sourcev29)
        event.getByLabel(pfShallowTag_labelv29, pfShallowTag_sourcev29)
        event.getByLabel(pfdeepbtag_labelv30, pfdeepbtag_sourcev30)
        event.getByLabel(pfShallowTag_labelv30, pfShallowTag_sourcev30)

        event.getByLabel(pfdeepbtag_labelv31, pfdeepbtag_sourcev31)
        event.getByLabel(pfShallowTag_labelv31, pfShallowTag_sourcev31)
        event.getByLabel(pfdeepbtag_labelv32, pfdeepbtag_sourcev32)
        event.getByLabel(pfShallowTag_labelv32, pfShallowTag_sourcev32)
        event.getByLabel(pfdeepbtag_labelv33, pfdeepbtag_sourcev33)
        event.getByLabel(pfShallowTag_labelv33, pfShallowTag_sourcev33)
        event.getByLabel(pfdeepbtag_labelv34, pfdeepbtag_sourcev34)
        event.getByLabel(pfShallowTag_labelv34, pfShallowTag_sourcev34)
        event.getByLabel(pfdeepbtag_labelv35, pfdeepbtag_sourcev35)
        event.getByLabel(pfShallowTag_labelv35, pfShallowTag_sourcev35)
        event.getByLabel(pfdeepbtag_labelv36, pfdeepbtag_sourcev36)
        event.getByLabel(pfShallowTag_labelv36, pfShallowTag_sourcev36)
        event.getByLabel(pfdeepbtag_labelv37, pfdeepbtag_sourcev37)
        event.getByLabel(pfShallowTag_labelv37, pfShallowTag_sourcev37)
        event.getByLabel(pfdeepbtag_labelv38, pfdeepbtag_sourcev38)
        event.getByLabel(pfShallowTag_labelv38, pfShallowTag_sourcev38)
        event.getByLabel(pfdeepbtag_labelv39, pfdeepbtag_sourcev39)
        event.getByLabel(pfShallowTag_labelv39, pfShallowTag_sourcev39)
        event.getByLabel(pfdeepbtag_labelv40, pfdeepbtag_sourcev40)
        event.getByLabel(pfShallowTag_labelv40, pfShallowTag_sourcev40)

        event.getByLabel(pfdeepbtag_labelv41, pfdeepbtag_sourcev41)
        event.getByLabel(pfShallowTag_labelv41, pfShallowTag_sourcev41)
        event.getByLabel(pfdeepbtag_labelv42, pfdeepbtag_sourcev42)
        event.getByLabel(pfShallowTag_labelv42, pfShallowTag_sourcev42)
        event.getByLabel(pfdeepbtag_labelv43, pfdeepbtag_sourcev43)
        event.getByLabel(pfShallowTag_labelv43, pfShallowTag_sourcev43)
        event.getByLabel(pfdeepbtag_labelv44, pfdeepbtag_sourcev44)
        event.getByLabel(pfShallowTag_labelv44, pfShallowTag_sourcev44)
        event.getByLabel(pfdeepbtag_labelv45, pfdeepbtag_sourcev45)
        event.getByLabel(pfShallowTag_labelv45, pfShallowTag_sourcev45)
        event.getByLabel(pfdeepbtag_labelv46, pfdeepbtag_sourcev46)
        event.getByLabel(pfShallowTag_labelv46, pfShallowTag_sourcev46)
        event.getByLabel(pfdeepbtag_labelv47, pfdeepbtag_sourcev47)
        event.getByLabel(pfShallowTag_labelv47, pfShallowTag_sourcev47)
        event.getByLabel(pfdeepbtag_labelv48, pfdeepbtag_sourcev48)
        event.getByLabel(pfShallowTag_labelv48, pfShallowTag_sourcev48)
        event.getByLabel(pfdeepbtag_labelv49, pfdeepbtag_sourcev49)
        event.getByLabel(pfShallowTag_labelv49, pfShallowTag_sourcev49)
        event.getByLabel(pfdeepbtag_labelv50, pfdeepbtag_sourcev50)
        event.getByLabel(pfShallowTag_labelv50, pfShallowTag_sourcev50)
        event.getByLabel(pfdeepbtag_labelv51, pfdeepbtag_sourcev51)
        event.getByLabel(pfShallowTag_labelv51, pfShallowTag_sourcev51)

        event.getByLabel(pfdeepbtag_labelv52, pfdeepbtag_sourcev52)
        event.getByLabel(pfShallowTag_labelv52, pfShallowTag_sourcev52)
        event.getByLabel(pfdeepbtag_labelv53, pfdeepbtag_sourcev53)
        event.getByLabel(pfShallowTag_labelv53, pfShallowTag_sourcev53)


        #print "getting offjets by label"
        event.getByLabel(offJets_label, offJets_source)
        event.getByLabel(offShallowTag_label, offShallowTag_source)
        #event.getByLabel(offJetsnoCHS_label, offJetsnoCHS_source)
        if runAOD:
            event.getByLabel(offbtag_label, offbtag_source)
            event.getByLabel(offdeepbtag_label, offdeepbtag_source)
            #event.getByLabel(offdeepbtag_bb_label, offdeepbtag_bb_source)
            #event.getByLabel(offdeepbtag_udsg_label, offdeepbtag_udsg_source)

        #Getting MET and HT handles
        #event.getByLabel(caloMet_label, caloMet_source)
        #event.getByLabel(caloMht_label, caloMht_source)
        #event.getByLabel(caloMhtNoPU_label, caloMhtNoPU_source)
        #event.getByLabel(pfMet_label, pfMet_source)
        #event.getByLabel(pfMht_label, pfMht_source)
        #event.getByLabel(offMet_label, offMet_source)

        #Getting Vertex handles
        event.getByLabel(FastPrimaryVertex_label, FastPrimaryVertex_source)
        event.getByLabel(FPVPixelVertices_label, FPVPixelVertices_source)
        event.getByLabel(PixelVertices_label, PixelVertices_source)
        event.getByLabel(VerticesPF_label, VerticesPF_source)
        event.getByLabel(VerticesL3_label, VerticesL3_source)
        event.getByLabel(VerticesOff_label, VerticesOff_source)

        #Getting Gen handles
        if isMC:
        #    event.getByLabel(genJets_label, genJets_source)
        #    if runAOD:
        #        event.getByLabel(genMet_label, genMet_source)
            event.getByLabel(genParticles_label, genParticles_source)
            event.getByLabel(generator_label, generator_source)
            event.getByLabel(pileUp_label, pileUp_source)


        event.getByLabel(triggerBitLabel, triggerBits)
        #print "finished getbylabel"

        #####################################################
        #####################################################
        
        names = event.object().triggerNames(triggerBits.product())
        triggerspassing = []
        for i,triggerName in enumerate(triggerNames):
            if triggerName.startswith("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_") or triggerName.startswith("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_"):
                index = names.triggerIndex(triggerName)
                if checkTriggerIndex(triggerName,index,names.triggerNames()):
                    triggerVars[triggerName][0] = triggerBits.product().accept(index)
                    #print triggerName,"acc:",triggerBits.product().accept(index)
                    #if triggerName.startswith("HLT") and not ( triggerName.startswith("NoFilter") or triggerName.endswith("FirstPath") or triggerName.endswith("FinalPath")):
                    if triggerName.startswith("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_") or triggerName.startswith("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_"):
                        if triggerBits.product().accept(index):
                            triggerspassing.append(triggerName)
                else:
                    triggerVars[triggerName][0] = 0

        # NOTE: Remove this if no trigger selection is required
        if doTriggerCut:
            if len(triggerspassing) < 1: 
                continue
        #print "finished triggers"
        ####################################################
        ####################################################

        run[0]          = event.eventAuxiliary().run()
        lumi[0]         = event.eventAuxiliary().luminosityBlock()
        evt[0]          = event.eventAuxiliary().event()

        if crun != run[0] or cls != lumi[0]:
            crun = run[0]
            cls = lumi[0]
            print "-------------- Processing: ",crun, cls," --------------"

        
        FastPrimaryVertex[0] = getVertex(FastPrimaryVertex_source)
        FPVPixelVertices[0] = getVertex(FPVPixelVertices_source)
        PixelVertices[0] = getVertex(PixelVertices_source)
        VerticesPF[0] = getVertex(VerticesPF_source)
        VerticesL3[0] = getVertex(VerticesL3_source)
        #VerticesOff[0]= getVertex(VerticesOff_source)
        vtx, nVtx = getVertices(VerticesOff_source)
        #print "saving"
        VerticesOff[0], nOffVertices[0] = vtx, nVtx 
        
        ####################################################
        ####################################################

        if isMC:
            trueVertex[0] = genParticles_source.productWithCheck().at(2).vertex().z()

        ####################################################
        ####################################################
        # Lepton Vectors
        offVertex = None
        if VerticesOff[0] > 0:
            offVertex = VerticesOff_source.productWithCheck().at(0)
            if not (offVertex.isFake() is False and offVertex.ndof() > 4 and abs(offVertex.z()) < 24 and abs(offVertex.position().Rho()) < 2):
                offVertex = None
                print "Offline Vertex did not pass the selection"

        FillElectronVector(offEle_source, offTightElectrons, idName )
        FillMuonVector(offMu_source, offTightMuons, offVertex, "loose")

        ####################################################
        ####################################################

        """
        #print "Filling MET and Stuff"
        FillVector(caloMet_source,caloMet, 0)
        FillVector(caloMht_source,caloMht, 0)
        FillVector(caloMhtNoPU_source,caloMhtNoPU, 0)
        FillVector(pfMet_source,pfMet, 0)
        FillVector(pfMht_source,pfMht, 0)

        l1Met[0],l1Met_phi[0],l1Mht[0],l1Mht_phi[0],l1HT[0] = -1,-1,-1,-1,-1
        for et in l1HT_source.productWithCheck():
            if et.getType()==ROOT.l1t.EtSum.kMissingEt:
                (l1Met[0],l1Met_phi[0]) = (et.et(),et.phi())
            elif et.getType()==ROOT.l1t.EtSum.kMissingHt:
                (l1Mht[0],l1Mht_phi[0]) = (et.et(),et.phi())
            elif et.getType()==ROOT.l1t.EtSum.kTotalEt:
                pass
            elif et.getType()==ROOT.l1t.EtSum.kTotalHt:
                l1HT[0] = et.et()

        FillVector(offMet_source,offMet, 0)

        if isMC:
            if runAOD:
                FillVector(genMet_source,genMet, 0)
        """
        ####################################################
        ##Edit this part to add the objects from the new paths to the tree
        ####################################################
        # Jets
        FillVector(caloJets_source,caloJets)
        FillVectorShallowTag(caloShallowTag_source, caloJets, jetVars, trkVars)


        FillVector(pfJets_source,pfJets)

        FillVectorShallowTag(pfShallowTag_source, pfJets, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv1)
        FillVectorShallowTag(pfShallowTag_sourcev1, pfJetsv1, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv2)
        FillVectorShallowTag(pfShallowTag_sourcev2, pfJetsv2, jetVars, trkVars)
        
        FillVector(pfJets_source,pfJetsv3)
        FillVectorShallowTag(pfShallowTag_sourcev3, pfJetsv3, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv4)
        FillVectorShallowTag(pfShallowTag_sourcev4, pfJetsv4, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv5)
        FillVectorShallowTag(pfShallowTag_sourcev5, pfJetsv5, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv6)
        FillVectorShallowTag(pfShallowTag_sourcev6, pfJetsv6, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv7)
        FillVectorShallowTag(pfShallowTag_sourcev7, pfJetsv7, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv8)
        FillVectorShallowTag(pfShallowTag_sourcev8, pfJetsv8, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv9)
        FillVectorShallowTag(pfShallowTag_sourcev9, pfJetsv9, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv10)
        FillVectorShallowTag(pfShallowTag_sourcev10, pfJetsv10, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv11)
        FillVectorShallowTag(pfShallowTag_sourcev11, pfJetsv11, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv12)
        FillVectorShallowTag(pfShallowTag_sourcev12, pfJetsv12, jetVars, trkVars)
        
        FillVector(pfJets_source,pfJetsv13)
        FillVectorShallowTag(pfShallowTag_sourcev13, pfJetsv13, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv14)
        FillVectorShallowTag(pfShallowTag_sourcev14, pfJetsv14, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv15)
        FillVectorShallowTag(pfShallowTag_sourcev15, pfJetsv15, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv16)
        FillVectorShallowTag(pfShallowTag_sourcev16, pfJetsv16, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv17)
        FillVectorShallowTag(pfShallowTag_sourcev17, pfJetsv17, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv18)
        FillVectorShallowTag(pfShallowTag_sourcev18, pfJetsv18, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv19)
        FillVectorShallowTag(pfShallowTag_sourcev19, pfJetsv19, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv20)
        FillVectorShallowTag(pfShallowTag_sourcev20, pfJetsv20, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv21)
        FillVectorShallowTag(pfShallowTag_sourcev21, pfJetsv21, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv22)
        FillVectorShallowTag(pfShallowTag_sourcev22, pfJetsv22, jetVars, trkVars)
        
        FillVector(pfJets_source,pfJetsv23)
        FillVectorShallowTag(pfShallowTag_sourcev23, pfJetsv23, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv24)
        FillVectorShallowTag(pfShallowTag_sourcev24, pfJetsv24, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv25)
        FillVectorShallowTag(pfShallowTag_sourcev25, pfJetsv25, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv26)
        FillVectorShallowTag(pfShallowTag_sourcev26, pfJetsv26, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv27)
        FillVectorShallowTag(pfShallowTag_sourcev27, pfJetsv27, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv28)
        FillVectorShallowTag(pfShallowTag_sourcev28, pfJetsv28, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv29)
        FillVectorShallowTag(pfShallowTag_sourcev29, pfJetsv29, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv30)
        FillVectorShallowTag(pfShallowTag_sourcev30, pfJetsv30, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv31)
        FillVectorShallowTag(pfShallowTag_sourcev31, pfJetsv31, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv32)
        FillVectorShallowTag(pfShallowTag_sourcev32, pfJetsv32, jetVars, trkVars)
        
        FillVector(pfJets_source,pfJetsv33)
        FillVectorShallowTag(pfShallowTag_sourcev33, pfJetsv33, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv34)
        FillVectorShallowTag(pfShallowTag_sourcev34, pfJetsv34, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv35)
        FillVectorShallowTag(pfShallowTag_sourcev35, pfJetsv35, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv36)
        FillVectorShallowTag(pfShallowTag_sourcev36, pfJetsv36, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv37)
        FillVectorShallowTag(pfShallowTag_sourcev37, pfJetsv37, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv38)
        FillVectorShallowTag(pfShallowTag_sourcev38, pfJetsv38, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv39)
        FillVectorShallowTag(pfShallowTag_sourcev39, pfJetsv39, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv40)
        FillVectorShallowTag(pfShallowTag_sourcev40, pfJetsv40, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv41)
        FillVectorShallowTag(pfShallowTag_sourcev41, pfJetsv41, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv42)
        FillVectorShallowTag(pfShallowTag_sourcev42, pfJetsv42, jetVars, trkVars)
        
        FillVector(pfJets_source,pfJetsv43)
        FillVectorShallowTag(pfShallowTag_sourcev43, pfJetsv43, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv44)
        FillVectorShallowTag(pfShallowTag_sourcev44, pfJetsv44, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv45)
        FillVectorShallowTag(pfShallowTag_sourcev45, pfJetsv45, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv46)
        FillVectorShallowTag(pfShallowTag_sourcev46, pfJetsv46, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv47)
        FillVectorShallowTag(pfShallowTag_sourcev47, pfJetsv47, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv48)
        FillVectorShallowTag(pfShallowTag_sourcev48, pfJetsv48, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv49)
        FillVectorShallowTag(pfShallowTag_sourcev49, pfJetsv49, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv50)
        FillVectorShallowTag(pfShallowTag_sourcev50, pfJetsv50, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv51)
        FillVectorShallowTag(pfShallowTag_sourcev51, pfJetsv51, jetVars, trkVars)

        FillVector(pfJets_source,pfJetsv52)
        FillVectorShallowTag(pfShallowTag_sourcev52, pfJetsv52, jetVars, trkVars)
        FillVector(pfJets_source,pfJetsv53)
        FillVectorShallowTag(pfShallowTag_sourcev53, pfJetsv53, jetVars, trkVars)







        #FillVector(l1Jets_source,l1Jets)
        FillVector(offJets_source,offJets,20, runAOD, offline = True, mc = isMC)
        FillVectorShallowTag(offShallowTag_source, offJets, jetVars, trkVars, debug=False)

        csvbtagTouples_calo = getBJetValuesforFilling(calobtag_source)
        deepcsvbtagTouples_calo = getBJetValuesforFilling(calodeepbtag_source)
        #FillBtag(calobtag_source, caloJets, caloJets.csv, caloJets.rankCSV,
        #         [CSVleadingCaloJet, CSVsecondCaloJet, CSVthirdCaloJet, CSVfourthCaloJet], nCSVCalogeZero, debug = False )
        FillBtagAlt(csvbtagTouples_calo, caloJets, caloJets.csv, caloJets.rankCSV, debug = False )
        #FillBtag(calodeepbtag_source, caloJets, caloJets.deepcsv, caloJets.rankDeepCSV,
        #         [DeepCSVleadingCaloJet, DeepCSVsecondCaloJet, DeepCSVthirdCaloJet, DeepCSVfourthCaloJet], nDeepCSVCalogeZero)
        FillBtagAlt(deepcsvbtagTouples_calo, caloJets, caloJets.deepcsv, caloJets.rankDeepCSV)
        ##FillBtag(calodeepbtag_bb_source, caloJets, caloJets.deepcsv_bb)
        ##FillBtag(calodeepbtag_udsg_source, caloJets, caloJets.deepcsv_udsg)
        #FillBtag(caloPUid_source, caloJets, caloJets.puId)

        #print "Filling pf btagging"
        csvbtagTouples_pf = getBJetValuesforFilling(pfbtag_source)
        deepcsvbtagTouples_pf = getBJetValuesforFilling(pfdeepbtag_source)
        deepcsvbtagTouples_pfv1 = getBJetValuesforFilling(pfdeepbtag_sourcev1)
        deepcsvbtagTouples_pfv2 = getBJetValuesforFilling(pfdeepbtag_sourcev2)
        deepcsvbtagTouples_pfv3 = getBJetValuesforFilling(pfdeepbtag_sourcev3)
        deepcsvbtagTouples_pfv4 = getBJetValuesforFilling(pfdeepbtag_sourcev4)
        deepcsvbtagTouples_pfv5 = getBJetValuesforFilling(pfdeepbtag_sourcev5)
        deepcsvbtagTouples_pfv6 = getBJetValuesforFilling(pfdeepbtag_sourcev6)
        deepcsvbtagTouples_pfv7 = getBJetValuesforFilling(pfdeepbtag_sourcev7)
        deepcsvbtagTouples_pfv8 = getBJetValuesforFilling(pfdeepbtag_sourcev8)
        deepcsvbtagTouples_pfv9 = getBJetValuesforFilling(pfdeepbtag_sourcev9)
        deepcsvbtagTouples_pfv10 = getBJetValuesforFilling(pfdeepbtag_sourcev10)

        deepcsvbtagTouples_pfv11 = getBJetValuesforFilling(pfdeepbtag_sourcev11)
        deepcsvbtagTouples_pfv12 = getBJetValuesforFilling(pfdeepbtag_sourcev12)
        deepcsvbtagTouples_pfv13 = getBJetValuesforFilling(pfdeepbtag_sourcev13)
        deepcsvbtagTouples_pfv14 = getBJetValuesforFilling(pfdeepbtag_sourcev14)
        deepcsvbtagTouples_pfv15 = getBJetValuesforFilling(pfdeepbtag_sourcev15)
        deepcsvbtagTouples_pfv16 = getBJetValuesforFilling(pfdeepbtag_sourcev16)
        deepcsvbtagTouples_pfv17 = getBJetValuesforFilling(pfdeepbtag_sourcev17)
        deepcsvbtagTouples_pfv18 = getBJetValuesforFilling(pfdeepbtag_sourcev18)
        deepcsvbtagTouples_pfv19 = getBJetValuesforFilling(pfdeepbtag_sourcev19)
        deepcsvbtagTouples_pfv20 = getBJetValuesforFilling(pfdeepbtag_sourcev20)

        deepcsvbtagTouples_pfv21 = getBJetValuesforFilling(pfdeepbtag_sourcev21)
        deepcsvbtagTouples_pfv22 = getBJetValuesforFilling(pfdeepbtag_sourcev22)
        deepcsvbtagTouples_pfv23 = getBJetValuesforFilling(pfdeepbtag_sourcev23)
        deepcsvbtagTouples_pfv24 = getBJetValuesforFilling(pfdeepbtag_sourcev24)
        deepcsvbtagTouples_pfv25 = getBJetValuesforFilling(pfdeepbtag_sourcev25)
        deepcsvbtagTouples_pfv26 = getBJetValuesforFilling(pfdeepbtag_sourcev26)
        deepcsvbtagTouples_pfv27 = getBJetValuesforFilling(pfdeepbtag_sourcev27)
        deepcsvbtagTouples_pfv28 = getBJetValuesforFilling(pfdeepbtag_sourcev28)
        deepcsvbtagTouples_pfv29 = getBJetValuesforFilling(pfdeepbtag_sourcev29)
        deepcsvbtagTouples_pfv30 = getBJetValuesforFilling(pfdeepbtag_sourcev30)

        deepcsvbtagTouples_pfv31 = getBJetValuesforFilling(pfdeepbtag_sourcev31)
        deepcsvbtagTouples_pfv32 = getBJetValuesforFilling(pfdeepbtag_sourcev32)
        deepcsvbtagTouples_pfv33 = getBJetValuesforFilling(pfdeepbtag_sourcev33)
        deepcsvbtagTouples_pfv34 = getBJetValuesforFilling(pfdeepbtag_sourcev34)
        deepcsvbtagTouples_pfv35 = getBJetValuesforFilling(pfdeepbtag_sourcev35)
        deepcsvbtagTouples_pfv36 = getBJetValuesforFilling(pfdeepbtag_sourcev36)
        deepcsvbtagTouples_pfv37 = getBJetValuesforFilling(pfdeepbtag_sourcev37)
        deepcsvbtagTouples_pfv38 = getBJetValuesforFilling(pfdeepbtag_sourcev38)
        deepcsvbtagTouples_pfv39 = getBJetValuesforFilling(pfdeepbtag_sourcev39)
        deepcsvbtagTouples_pfv40 = getBJetValuesforFilling(pfdeepbtag_sourcev40)

        deepcsvbtagTouples_pfv41 = getBJetValuesforFilling(pfdeepbtag_sourcev41)
        deepcsvbtagTouples_pfv42 = getBJetValuesforFilling(pfdeepbtag_sourcev42)
        deepcsvbtagTouples_pfv43 = getBJetValuesforFilling(pfdeepbtag_sourcev43)
        deepcsvbtagTouples_pfv44 = getBJetValuesforFilling(pfdeepbtag_sourcev44)
        deepcsvbtagTouples_pfv45 = getBJetValuesforFilling(pfdeepbtag_sourcev45)
        deepcsvbtagTouples_pfv46 = getBJetValuesforFilling(pfdeepbtag_sourcev46)
        deepcsvbtagTouples_pfv47 = getBJetValuesforFilling(pfdeepbtag_sourcev47)
        deepcsvbtagTouples_pfv48 = getBJetValuesforFilling(pfdeepbtag_sourcev48)
        deepcsvbtagTouples_pfv49 = getBJetValuesforFilling(pfdeepbtag_sourcev49)
        deepcsvbtagTouples_pfv50 = getBJetValuesforFilling(pfdeepbtag_sourcev50)
        deepcsvbtagTouples_pfv51 = getBJetValuesforFilling(pfdeepbtag_sourcev51)

        deepcsvbtagTouples_pfv52 = getBJetValuesforFilling(pfdeepbtag_sourcev52)
        deepcsvbtagTouples_pfv53 = getBJetValuesforFilling(pfdeepbtag_sourcev53)



        #FillBtag(pfbtag_source, pfJets, pfJets.csv, pfJets.rankCSV,
        #         [CSVleadingPFJet, CSVsecondPFJet, CSVthirdPFJet, CSVfourthPFJet], nCSVPFgeZero)
        FillBtagAlt(csvbtagTouples_pf, pfJets, pfJets.csv, pfJets.rankCSV)
        #FillBtag(pfdeepbtag_source, pfJets, pfJets.deepcsv, pfJets.rankDeepCSV,
        #         [DeepCSVleadingPFJet, DeepCSVsecondPFJet, DeepCSVthirdPFJet, DeepCSVfourthPFJet], nDeepCSVPFgeZero, debug = False)
        FillBtagAlt(deepcsvbtagTouples_pf, pfJets, pfJets.deepcsv, pfJets.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv1, pfJetsv1, pfJetsv1.deepcsv, pfJetsv1.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv2, pfJetsv2, pfJetsv2.deepcsv, pfJetsv2.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv3, pfJetsv3, pfJetsv3.deepcsv, pfJetsv3.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv4, pfJetsv4, pfJetsv4.deepcsv, pfJetsv4.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv5, pfJetsv5, pfJetsv5.deepcsv, pfJetsv5.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv6, pfJetsv6, pfJetsv6.deepcsv, pfJetsv6.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv7, pfJetsv7, pfJetsv7.deepcsv, pfJetsv7.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv8, pfJetsv8, pfJetsv8.deepcsv, pfJetsv8.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv9, pfJetsv9, pfJetsv9.deepcsv, pfJetsv9.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv10, pfJetsv10, pfJetsv10.deepcsv, pfJetsv10.rankDeepCSV)

        FillBtagAlt(deepcsvbtagTouples_pfv11, pfJetsv11, pfJetsv11.deepcsv, pfJetsv11.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv12, pfJetsv12, pfJetsv12.deepcsv, pfJetsv12.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv13, pfJetsv13, pfJetsv13.deepcsv, pfJetsv13.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv14, pfJetsv14, pfJetsv14.deepcsv, pfJetsv14.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv15, pfJetsv15, pfJetsv15.deepcsv, pfJetsv15.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv16, pfJetsv16, pfJetsv16.deepcsv, pfJetsv16.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv17, pfJetsv17, pfJetsv17.deepcsv, pfJetsv17.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv18, pfJetsv18, pfJetsv18.deepcsv, pfJetsv18.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv19, pfJetsv19, pfJetsv19.deepcsv, pfJetsv19.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv20, pfJetsv20, pfJetsv20.deepcsv, pfJetsv20.rankDeepCSV)

        FillBtagAlt(deepcsvbtagTouples_pfv21, pfJetsv21, pfJetsv21.deepcsv, pfJetsv21.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv22, pfJetsv22, pfJetsv22.deepcsv, pfJetsv22.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv23, pfJetsv23, pfJetsv23.deepcsv, pfJetsv23.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv24, pfJetsv24, pfJetsv24.deepcsv, pfJetsv24.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv25, pfJetsv25, pfJetsv25.deepcsv, pfJetsv25.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv26, pfJetsv26, pfJetsv26.deepcsv, pfJetsv26.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv27, pfJetsv27, pfJetsv27.deepcsv, pfJetsv27.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv28, pfJetsv28, pfJetsv28.deepcsv, pfJetsv28.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv29, pfJetsv29, pfJetsv29.deepcsv, pfJetsv29.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv30, pfJetsv30, pfJetsv30.deepcsv, pfJetsv30.rankDeepCSV)

        FillBtagAlt(deepcsvbtagTouples_pfv31, pfJetsv31, pfJetsv31.deepcsv, pfJetsv31.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv32, pfJetsv32, pfJetsv32.deepcsv, pfJetsv32.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv33, pfJetsv33, pfJetsv33.deepcsv, pfJetsv33.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv34, pfJetsv34, pfJetsv34.deepcsv, pfJetsv34.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv35, pfJetsv35, pfJetsv35.deepcsv, pfJetsv35.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv36, pfJetsv36, pfJetsv36.deepcsv, pfJetsv36.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv37, pfJetsv37, pfJetsv37.deepcsv, pfJetsv37.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv38, pfJetsv38, pfJetsv38.deepcsv, pfJetsv38.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv39, pfJetsv39, pfJetsv39.deepcsv, pfJetsv39.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv40, pfJetsv40, pfJetsv40.deepcsv, pfJetsv40.rankDeepCSV)

        FillBtagAlt(deepcsvbtagTouples_pfv41, pfJetsv41, pfJetsv41.deepcsv, pfJetsv41.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv42, pfJetsv42, pfJetsv42.deepcsv, pfJetsv42.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv43, pfJetsv43, pfJetsv43.deepcsv, pfJetsv43.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv44, pfJetsv44, pfJetsv44.deepcsv, pfJetsv44.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv45, pfJetsv45, pfJetsv45.deepcsv, pfJetsv45.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv46, pfJetsv46, pfJetsv46.deepcsv, pfJetsv46.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv47, pfJetsv47, pfJetsv47.deepcsv, pfJetsv47.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv48, pfJetsv48, pfJetsv48.deepcsv, pfJetsv48.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv49, pfJetsv49, pfJetsv49.deepcsv, pfJetsv49.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv50, pfJetsv50, pfJetsv50.deepcsv, pfJetsv50.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv51, pfJetsv51, pfJetsv51.deepcsv, pfJetsv51.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv52, pfJetsv52, pfJetsv52.deepcsv, pfJetsv52.rankDeepCSV)
        FillBtagAlt(deepcsvbtagTouples_pfv53, pfJetsv53, pfJetsv53.deepcsv, pfJetsv53.rankDeepCSV)


        #FillBtag(pfdeepbtag_bb_source, pfJets, pfJets.deepcsv_bb)
        #FillBtag(pfdeepbtag_udsg_source, pfJets, pfJets.deepcsv_udsg)
        

        if runAOD:
            FillBtag(offbtag_source, offJets, offJets.csv, offJets.rankCSV,
                     [CSVleadingOffJet, CSVsecondOffJet, CSVthirdOffJet, CSVfourthOffJet])#, nCSVOffgeZero)
            FillBtag(offdeepbtag_source, offJets, offJets.deepcsv_b)#, nDeepCSVOffgeZero)
            FillBtag(offdeepbtag_bb_source, offJets, offJets.deepcsv_bb)
            FillBtag(offdeepbtag_udsg_source, offJets, offJets.deepcsv_udsg)
        else:
            makeCSVRanking(offJets, offJets.csv, offJets.rankCSV)

            
        makeDeepCSVSumRanking(offJets, offJets.deepcsv, offJets.deepcsv_b, offJets.deepcsv_bb, offJets.rankDeepCSV)
        LeptonOverlap(offJets, offTightMuons, offTightElectrons, offJets.lepOverlap04Tight)

        """ Removed GenJets for now 
        if isMC:
            FillVector(genJets_source,genJets,15)
        """

        resetArray(offJets.matchCalo, -1)
        resetArray(offJets.matchPF, -1)
        for i in range(caloJets.num[0]):
            caloJets.matchOff[i], dR = Matching(caloJets.phi[i],caloJets.eta[i],offJets)
            offJets.matchCalo[int(caloJets.matchOff[i])] = i
            caloJets.matchGen[i] = -1

        #print "Matching pf jets to off jets"
        for i in range(pfJets.num[0]):
            pfJets.matchOff[i], dR = Matching(pfJets.phi[i],pfJets.eta[i],offJets)
            offJets.matchPF[int(pfJets.matchOff[i])] = i
            pfJets.matchGen[i] = -1
            #print "Online PF jet (pt = ",pfJets.pt[i],") ",i," matched to offline jet (pt = ",offJets.pt[pfJets.matchOff[i]],") ",pfJets.matchOff[i]," with dR = ",dR
            
        #print "Matching l1 jets to off jets"
        #for i in range(l1Jets.num[0]):
        #    l1Jets.matchOff[i], dR = Matching(l1Jets.phi[i],l1Jets.eta[i],offJets)
        #    l1Jets.matchGen[i] = -1

        """ Removed GenJets for now 
        if isMC:
            ####################################################
            # Gen patricle for Jet
            for i in range(genJets.num[0]):
                genJets.mcFlavour[i] = -100
                genJets.mcPt[i] = -100

            for i in range(caloJets.num[0]):
                caloJets.matchGen[i], dR = Matching(caloJets.phi[i],caloJets.eta[i],genJets)

            #Matching pf jets to off and gen jets
            for i in range(pfJets.num[0]):
                pfJets.matchGen[i], dR = Matching(pfJets.phi[i],pfJets.eta[i],genJets)
                
                
            #Matching l1 jets to off and gen jets
            #for i in range(l1Jets.num[0]):
            #    l1Jets.matchGen[i], dR = Matching(l1Jets.phi[i],l1Jets.eta[i],genJets)
        
            #Matching gen jets to off jets
            #for i in range(offJets.num[0]):
            #    offJets.matchGen[i], dR = Matching(offJets.phi[i],offJets.eta[i],genJets)

            for genParticle in genParticles_source.productWithCheck():
                if genParticle.pt()<5: continue
                if not (abs(genParticle.pdgId()) in [21,1,2,3,4,5,11,13,15]): continue
                if genParticle.mother().pt()>5 and (abs(genParticle.mother().pdgId()) in [21,1,2,3,4,5,11,13]): continue
                if evt[0]==7826939:
                    print "genParticle:"
                    print genParticle.pt(),genParticle.eta(),genParticle.phi(),genParticle.pdgId()
                    print "genJets:"
                for i in range(genJets.num[0]):
                    if genParticle.pt()<0.2*genJets.pt[i]: continue
                    if deltaR(genParticle.eta(),genParticle.phi(),genJets.eta[i],genJets.phi[i])<0.4:
                        if evt[0]==7826939:
                            print genJets.pt[i],genJets.eta[i],genJets.eta[i],genJets.mcFlavour[i]
                            print "not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1]):",not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1])
                        if abs(genParticle.pdgId())==5:
                            if genJets.mcFlavour[i]!=5 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==4 and not int(abs(genJets.mcFlavour[i])) in [5]:
                            if genJets.mcFlavour[i]!=4 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==3 and not int(abs(genJets.mcFlavour[i])) in [5,4]:
                            if genJets.mcFlavour[i]!=3 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==2 and not int(abs(genJets.mcFlavour[i])) in [5,4,3]:
                            if genJets.mcFlavour[i]!=2 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==1 and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2]:
                            if genJets.mcFlavour[i]!=1 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId())==21 and not (int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1]):
                            if genJets.mcFlavour[i]!=21 or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [11,13] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21]:
                            if not (genJets.mcFlavour[i] in [11,13]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [15] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21,11,13]:
                            if not (genJets.mcFlavour[i] in [15]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        elif abs(genParticle.pdgId()) in [22] and not int(abs(genJets.mcFlavour[i])) in [5,4,3,2,1,21,11,13,15]:
                            if not (genJets.mcFlavour[i] in [22]) or genParticle.pt()>genJets.mcPt[i]:
                                genJets.mcFlavour[i] = genParticle.pdgId()
                                genJets.mcPt[i]      = genParticle.pt()
                        if evt[0]==7826939:
                            print "newFlav:",genJets.mcFlavour[i]
            
            FillMCFlavour(offJets, offJets.matchGen, genJets.mcFlavour, offJets.mcFlavour)
            #FillMCFlavour(offCSVJets, offCSVJets.matchGen, genJets.mcFlavour, offCSVJets.mcFlavour)
            #FillMCFlavour(offDeepCSVJets, offDeepCSVJets.matchGen, genJets.mcFlavour, offDeepCSVJets.mcFlavour)
            FillMCFlavour(caloJets, caloJets.matchGen, genJets.mcFlavour, caloJets.mcFlavour)
            FillMCFlavour(pfJets, pfJets.matchGen, genJets.mcFlavour, pfJets.mcFlavour)
        """
        ####################################################
        ####################################################
        #print "Making cleanCollections"

        #cleanCollection(offJets, offCleanJets, "lepOverlap04Tight", 0.0, offCleanJets.offOrder, verbose = False)
        #cleanCollection(offJets, offClean05Jets, "lepOverlap05Tight", 0.0, offClean05Jets.offOrder)

        #sortJetCollection(offCleanJets, offCleanCSVJets, "csv", offCleanCSVJets.rankpt)
        #sortJetCollection(offCleanJets, offCleanDeepCSVJets, "deepcsv", offCleanDeepCSVJets.rankpt)
        
        if isMC:
            if bunchCrossing>=pileUp_source.productWithCheck().size() or pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing()!=0:
                found=False
                for bunchCrossing in range(pileUp_source.productWithCheck().size()):
                    if pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing() == 0 :
                        found=True;
                        break
                if not found:
                    Exception("Check pileupSummaryInfos!")
                print "I'm using bunchCrossing=",bunchCrossing
            pu[0] = pileUp_source.productWithCheck().at(bunchCrossing).getTrueNumInteractions()
            print pu[0], pileUp_source.productWithCheck().at(bunchCrossing).getTrueNumInteractions()
            """
            wPURunC[0] = getPUweight("RunC", pu[0])
            wPURunD[0] = getPUweight("RunD", pu[0])
            wPURunE[0] = getPUweight("RunE", pu[0])
            wPURunF[0] = getPUweight("RunF", pu[0])
            wPURunCF[0] = getPUweight("RunC-F", pu[0])
            """
            ptHat[0]    = generator_source.product().qScale()

            maxPUptHat[0] = -1
            for ptHat_ in pileUp_source.productWithCheck().at(bunchCrossing).getPU_pT_hats():
                maxPUptHat[0] = max(maxPUptHat[0],ptHat_)
                
        if iev%10000==1: print "Event: ",iev," done."
        #print "Event: ",iev," done."
        nPassHisto.Fill(1)
        tree.Fill()
        
    f.Write()
    f.Close()
    dir_ = os.getcwd()
    print "Total time: {0:10f} s".format(time.time()-t0)
    print "Filesize of {0:8f} MB".format(os.path.getsize(dir_+"/"+fileOutput) * 1e-6)

if __name__ == "__main__":
    secondaryFiles = [
        #"file:/afs/cern.ch/work/k/koschwei/public/TTTo2L2Nu_RunIIFall17DRPremix_94X_TSG_GEN-SIM-RAW_EC33C288-051E-E811-95B7-1866DA7F98A8.root"
        'file:/mnt/t3nfs01/data01/shome/koschwei/scratch/MuonEGRunC_RAW_300107_348E3CF3-6974-E711-80DE-02163E01A5DC.root'
        #"/store/mc/RunIIFall17DRPremix/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RAW/TSG_94X_mc2017_realistic_v11-v1/30000/14E63E79-041E-E811-A72A-3417EBE5289A.root"
        #"file:/mnt/t3nfs01/data01/shome/koschwei/scratch/MuonEGRunC_RAW_300107_348E3CF3-6974-E711-80DE-02163E01A5DC.root"
        #"/store/data/Run2017C/MuonEG/RAW/v1/000/300/155/00000/5219D654-9176-E711-99E5-02163E019DED.root"
    ]
    filesInput = [
        #"file:/afs/cern.ch/work/k/koschwei/public/TTTo2L2Nu_RunIIFall17MiniAOD_94X_TSG_MINIAOD_A25F33EE-5A1F-E811-ABCC-90B11C48DA2F.root"
        #"file:/mnt/t3nfs01/data01/shome/koschwei/scratch/MuonEGRunC_MiniAODv2_ReReco_300107_4AABB4B6-5037-E811-9F4A-0CC47A5FBE35.root"
        #"/store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/TSG_94X_mc2017_realistic_v11-v1/30000/AEA25818-FB21-E811-BBFB-848F69FD444B.root"
        "file:/afs/cern.ch/work/k/koschwei/public/MuonEGRunC_MiniAOD_ReReco_300107_221A60AF-A1EC-E711-9510-A4BF011259A8.root"
        #"file:/mnt/t3nfs01/data01/shome/koschwei/scratch/MuonEGRunC_MiniAODv2_ReReco_300107_4AABB4B6-5037-E811-9F4A-0CC47A5FBE35.root"
        #"/store/data/Run2017C/MuonEG/MINIAOD/PromptReco-v2/000/300/155/00000/14B10372-AF77-E711-A0A0-02163E0133CC.root"
    ]
    fileOutput = "tree_phase1.root"
    maxEvents = 10
    launchNtupleFromHLT(fileOutput,filesInput,secondaryFiles,maxEvents, local = True, preProcessing=False)

    
