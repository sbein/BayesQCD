#Welcome to the industrial age of Sam's rebalance and smear code. You're going to have a lot of fun!
import os,sys
from ROOT import *
from array import array
from glob import glob
from utils import *
#from ra2blibs import *
import time

###stuff that would be nice in a config file
mhtjetetacut = 5.0 # also needs be be changed in UsefulJet.h
AnHardMetJetPtCut = 30.0
cutoff = 15.0
isdata = False
rebalancedMetCut = 90
#rebalancedMetCut = 150


debugmode = False
searchmode = False
usefullprior = False

##load in delphes libraries to access input collections:
gSystem.Load("/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/SampleProduction/delphes/build/libDelphes")
##load in UsefulJet class, which the rebalance and smear code uses
gROOT.ProcessLine(open('src/UsefulJet.cc').read())
exec('from ROOT import *')
gROOT.ProcessLine(open('src/BayesRandS.cc').read())
exec('from ROOT import *')


##read in command line arguments
defaultInfile_ = "../SampleProduction/delphes/output_ak5/delphes_gjet7.root"
#T2qqGG.root
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", type=int, default=1,help="analyzer script to batch")
parser.add_argument("-printevery", "--printevery", type=int, default=100,help="short run")
parser.add_argument("-fin", "--fnamekeyword", type=str,default=defaultInfile_,help="file")
parser.add_argument("-bootstrap", "--bootstrap", type=str, default='0',help="boot strapping (0,1of5,2of5,3of5,...)")
parser.add_argument("-quickrun", "--quickrun", type=bool, default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
bootstrap = args.bootstrap
verbosity = args.verbosity
printevery = args.printevery
quickrun = args.quickrun
if quickrun: 
    n2process = 50000
    if 'T2' in fnamekeyword: n2process = 5000
else: n2process = 9999999999999
mktree = False
BTag_Cut = 0.5 #Delphes b-tagging binary
if bootstrap=='0': 
    bootstrapmode = False
    bootupfactor = 1
else: 
    bootstrapmode = True
    from random import randint
    thisbootstrap, nbootstraps = bootstrap.split('of')
    thisbootstrap, nbootstraps = int(thisbootstrap), int(nbootstraps)
    print 'thisbootstrap, nbootstraps', thisbootstrap, nbootstraps
    bootupfactor = nbootstraps


#Dictionary list of signal regions
regionCuts = {}
pi = 3.14159
Inf = 9999
#varlist =                    ['St',    'Ht',     'HardMet','NJets','BTags','NPhotons', 'PhoPt',  'PhoEta','MinDeltaPhi','St5OverSt']
#regionCuts['NoCuts']        = [[0,Inf], [0,Inf], [0,Inf], [0,Inf],[0,Inf],  [0,Inf],  [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,1.03]]
regionCuts['Baseline0Pho']   = [[300,Inf],[300,Inf],[300,Inf],[2,Inf],[0,Inf],  [0,0],    [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]]



inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/pythia_based/T1bbbb-1800-200_delphesout.root']
inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/madgraph_based/archive/T1bbbb_1800_200_B_delphesout.root']
inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/pythia_based/T1qqqq-1300-100_delphesout.root']#pythia
#inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/pythia_based']
#inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/pythia_based/T1bbbb-1300-1100_delphesout.root']
inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/madgraph_based/T1bbbb_1300_1100_delphesout.root']
#inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/madgraph_based/archive/T5qqqqVV_1800_100_delphesout.root']
inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/madgraph_based/T5qqqqVV_1800_100_delphesout.root']
inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/madgraph_based/T5qqqqVV_1800_100_delphesout.root']
#inputFiles = ['/nfs/dust/cms/user/mrowietm/delphesfiles/pythia_based/T5qqqqVV-1800-100_delphesout.root']


##declare and load a tree
c = TChain('Delphes')
for fname in inputFiles: 
    print 'adding', fname
    c.Add(fname)
nentries = c.GetEntries()
c.Show(0)
n2process = min(n2process, nentries)
print 'n(entries) =', n2process

##feed the tree to delphes, set up which branches need to be used
treeReader = ExRootTreeReader(c)
numberOfEntries = treeReader.GetEntries()
branchHT = treeReader.UseBranch("ScalarHT")
branchJets = treeReader.UseBranch("Jet")
branchGenJet = treeReader.UseBranch("GenJet")
branchMissingET = treeReader.UseBranch("MissingET")
branchGenMissingET = treeReader.UseBranch("GenMissingET")
branchPhoton = treeReader.UseBranch("Photon")
branchElectron = treeReader.UseBranch("Electron")
branchEvent = treeReader.UseBranch("Event")
branchMuon = treeReader.UseBranch("Muon")

#branchTower = treeReader.UseBranch("Tower")
#branchEFlowTower = treeReader.UseBranch("EFlowPhoton")


branchParticles = treeReader.UseBranch("Particle")


varlist = ['St', 'Ht', 'HardMet','NJets','BTags','NPhotons', 'PhoPt', 'PhoEta','MinDeltaPhi','St5OverSt']
indexVar = {}
for ivar, var in enumerate(varlist): indexVar[var] = ivar
indexVar[''] = -1
nmain = len(varlist)

def selectionFeatureVector(fvector, regionkey='', omitcuts='', omitcuts_dphi=''):
    iomits, iomits_dphi = [], []  
    for cut in omitcuts.split('Vs'): iomits.append(indexVar[cut])
    for i, feature in enumerate(fvector):
        if i==nmain: break
        if i in iomits: continue
        if not (feature>=regionCuts[regionkey][i][0] and feature<=regionCuts[regionkey][i][1]): 
            return False
    return True


templateFileName = 'usefulthings/llhd-prior-coarse-ak4.root'
templateFileName = 'usefulthings/llhd-prior-coarse-ak4-extra.root'

templateFileName = 'usefulthings/llhd-prior-coarse-doublewidth.root'##need to re-run with new templates
templateFileName = 'usefulthings/llhd-prior-coarse-240files.root'
templateFileName = 'usefulthings/llhd-prior-coarse-1p2width.root'


#templateFileName = 'usefulthings/llhd-prior-MetSignificance.root'
ftemplate = TFile(templateFileName)
print 'using templates from',templateFileName
hPtTemplate = ftemplate.Get('hPtTemplate')
templatePtAxis = hPtTemplate.GetXaxis()
hEtaTemplate = ftemplate.Get('hEtaTemplate')
templateEtaAxis = hEtaTemplate.GetXaxis()



#option for using FullSim-based prior
if usefullprior:
    priorFileName = templateFileName
    priorFileName = 'usefulthings/ResponseFunctionsMC17AllFilters_deepCsv.root'
    fprior = TFile(priorFileName)
else:
    fprior = TFile(templateFileName)

hHtTemplate = fprior.Get('hHtTemplate')
templateStAxis = hHtTemplate.GetXaxis()


##Create output file
infileID = fnamekeyword.split('/')[-1].replace('.root','')
newname = 'ma5plots.root'
fnew = TFile(newname, 'recreate')
print 'creating', newname
if mktree:
    treefile = TFile('littletreeLowHardMet'+fnamekeyword+'.root','recreate')
    littletree = TTree('littletree','littletree')
    prepareLittleTree(littletree)

hSt = TH1F('hSt','hSt',120,0,2500)
hSt.Sumw2()
hStWeighted = TH1F('hStWeighted','hStWeighted',120,0,2500)
hStWeighted.Sumw2()
hGenMetGenHardMetRatio = TH1F('hGenMetGenHardMetRatio','hGenMetGenHardMetRatio',50,0,5)
hGenMetGenHardMetRatio.Sumw2()
hPassFit = TH1F('hPassFit','hPassFit',10,0,10)
hPassFit.Sumw2()
hTotFit = TH1F('hTotFit','hTotFit',10,0,10)
hTotFit.Sumw2()



hPtGenLeptons = TH1F('hPtGenLeptons','hPtGenLeptons',25,0,500)
histoStyler(hPtGenLeptons,kBlack)
hEtaGenLeptons = TH1F('hEtaGenLeptons','hEtaGenLeptons',50,-5,5)
histoStyler(hEtaGenLeptons,kBlack)

hPtRecoLeptons = TH1F('hPtRecoLeptons','hPtRecoLeptons',25,0,500)
histoStyler(hPtRecoLeptons,kBlack)
hEtaRecoLeptons = TH1F('hEtaRecoLeptons','hEtaRecoLeptons',50,-5,5)
histoStyler(hEtaRecoLeptons,kBlack)

hPtRecoIsoLeptons = TH1F('hPtRecoIsoLeptons','hPtRecoIsoLeptons',25,0,500)
histoStyler(hPtRecoIsoLeptons,kBlack)
hEtaRecoIsoLeptons = TH1F('hEtaRecoIsoLeptons','hEtaRecoIsoLeptons',50,-5,5)
histoStyler(hEtaRecoIsoLeptons,kBlack)

hNLepton = TH1F('hNLepton','hNLepton',5,0,5)
histoStyler(hNLepton,kBlack)

#GleanTemplatesFromFile(ftemplate)
print 'listing contents of prior:'
fprior.ls('splines/*')

#GleanTemplatesFromFile(ftemplate, fprior)
GleanTemplatesFromFile(fprior)


binning['Ht']=[30,0,3000]

histoStructDict = {}
for region in regionCuts:
    for var in varlist:
        histname = region+'_'+var
        histoStructDict[histname] = mkHistoStruct(histname, binning)

xsec_times_lumi_over_nevents = 1.0
t0 = time.time()


#n2process = 10000

for ientry in range(n2process):


    if ientry%printevery==0:
        print "processing event", ientry, '/', n2process, 'time', time.time()-t0

    if debugmode:
        #if not ientry>122000: continue
        if not ientry in [5068]: continue
        print 'top of event', ientry, '='*20
        a = 2
    c.GetEntry(ientry)


    weight = xsec_times_lumi_over_nevents

    #if not branchMissingET[0].MET>175: continue

    acme_objects = vector('TLorentzVector')()
    recophotons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using



    #######if not nspikes==0: continue #####spikes!
    ###if not len(recophotons)==0: continue

    if len(recophotons)>0: 
        phopt, phoeta = recophotons[0].Pt(), recophotons[0].Eta()
    else:
        phopt, phoeta = -1, -11

    genweight = branchEvent[0].Weight
    if debugmode or searchmode:
        genthingies = []
        rawgenmhtthingy = TLorentzVector()
        for igp_, gp_ in enumerate(branchParticles):
            if not gp_.Status==1: continue
            if abs(gp_.PID) in [12,14,16,1000022]: continue
            thistlv = TLorentzVector()
            thistlv.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
            rawgenmhtthingy-=thistlv
            genthingies.append([thistlv, igp_, gp_.PID])

            if abs(gp_.PID) in [11,13] and gp_.PT>5: 
                #print 'filling electron or muon', gp_.PID
                fillth1(hPtGenLeptons,gp_.PT)
                fillth1(hEtaGenLeptons,gp_.Eta)                
            if not gp_.PT>30: continue

            #print 'weve gota a gen', gp_.PID, 'pT =', gp_.PT, 'Eta =', gp_.Eta, gp_.Phi

    if len(recophotons)>1: dphi = abs(recophotons[0].DeltaPhi(recophotons[1]))
    else: dphi = -1

    #print 'now processing a', ientry

    recoelectrons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for iel, el in enumerate(branchElectron):
        if not el.PT>10: continue
        fillth1(hPtRecoLeptons, el.PT)        
        if not el.IsolationVar<0.2: continue
        fillth1(hPtRecoIsoLeptons, el.PT)        
        tlvel = TLorentzVector()
        tlvel.SetPtEtaPhiE(el.PT, el.Eta, el.Phi, el.PT*TMath.CosH(el.Eta))
        if debugmode:
            print ientry, 'acme electron pt and iso', el.PT, el.IsolationVar
        if el.PT>10: acme_objects.push_back(tlvel)		
        if not abs(el.Eta)<2.4: continue		
        recoelectrons.push_back(tlvel)


    recomuons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for imu, mu in enumerate(branchMuon):
        if not mu.PT>10: continue
        fillth1(hPtRecoLeptons, mu.PT)
        if not mu.IsolationVar<0.2: continue
        fillth1(hPtRecoIsoLeptons, mu.PT)
        tlvmu = TLorentzVector()
        tlvmu.SetPtEtaPhiE(mu.PT, mu.Eta, mu.Phi, mu.PT*TMath.CosH(mu.Eta))
        if debugmode:
            print ientry, 'acme muon pt and iso', mu.PT, mu.IsolationVar			
        if mu.PT>10: acme_objects.push_back(tlvmu)			
        if not abs(mu.Eta)<5.0: continue		
        recomuons.push_back(tlvmu)		

    fillth1(hNLepton, len(recomuons)+len(recoelectrons))
    if not len(recoelectrons)==0: continue    
    if not len(recomuons)==0: continue

    AcmeVector = TLorentzVector()
    AcmeVector.SetPxPyPzE(0,0,0,0)
    for obj in acme_objects: AcmeVector+=obj		

    _Templates_.AcmeVector = AcmeVector	

    ##declare empty vector of UsefulJets (in c++, std::vector<UsefulJet>):
    recojets = vector('UsefulJet')()


    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using

    onlygoodjets = True
    for ijet, jet in enumerate(branchJets):
    ##for ijet, jet in enumerate(branchGenJet):        
    #if debugmode: print 'initially looking at jet', ijet, jet.PT, jet.EhadOverEem
        if not jet.PT>15: continue
        if not abs(jet.Eta)<5: continue

        
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))

        if debugmode:
            if abs(jet.PT-839.98)<1: 
                gmethods = dir(jet)
                print gmethods
                for ipart, part in enumerate(jet.Particles):
                    print ijet, ipart, 'jetparticle', ipart, part.PT, part.PID, part.Eta
                part0 = jet.Particles[0]
                part1 = jet.Particles[1]
                tlv0 = TLorentzVector()
                tlv1 = TLorentzVector()
                
                tlv0.SetPtEtaPhiE(part0.PT, part0.Eta, part0.Phi, part0.E)
                tlv1.SetPtEtaPhiE(part1.PT, part1.Eta, part1.Phi, part1.E)                
                print 'deltaR', tlv0.DeltaR(tlv1)
                for gmethod in gmethods[:40]:
                    try:  print 'method', gmethod+':', getattr(jet, gmethod)
                    except: pass  
                for igjet, gjet in enumerate(branchGenJet):		
                    #if not jet.PT>15: continue
                    #if not abs(jet.Eta)<5: continue
                    tlvgjet = TLorentzVector()
                    tlvgjet.SetPtEtaPhiE(gjet.PT, gjet.Eta, gjet.Phi, gjet.PT*TMath.CosH(gjet.Eta))
                    print 'testing this gen jet with pT', gjet.PT, tlvjet.DeltaR(tlvgjet)
                        
        ##ujet = UsefulJet(tlvjet, jet.BTag)
        ujet = UsefulJet(tlvjet, 0)
        if calcMinDr(acme_objects, ujet, 0.3)<0.3:
            if debugmode: print 'throwing away reco jet', ujet.Pt()
            continue
        if calcMinDr(acme_objects, ujet, 0.3)<0.55:
            onlygoodjets = False

        recojets.push_back(ujet)
        ##if ujet.Pt()>50 and jet.EhadOverEem>100: 
        ##    onlygoodjets = False
        ##    #print ijet, 'bad jet'
    if not onlygoodjets: 
        a = 2
        #print 'skipping on jet requirement'
        #continue


    tHardMetVec = getHardMet(recojets,AnHardMetJetPtCut, mhtjetetacut) # for debugging
    tHardMhtVec = tHardMetVec.Clone() # for debugging
    tHardMetVec-=AcmeVector # for debugging
    tHardMetPt, tHardMetPhi = tHardMetVec.Pt(), tHardMetVec.Phi() # for debugging
    tHardMhtPt, tHardMhtPhi = tHardMhtVec.Pt(), tHardMhtVec.Phi() # for debugging

    #print 'now processing b', ientry

    if debugmode:
        print 'len(acme_objects)', len(acme_objects)

    genmet = TLorentzVector()
    allgenjets = vector('UsefulJet')()
    genjets_ = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #user should interface with the input format of the event using
    for ijet, jet in enumerate(branchGenJet):		
        #if not jet.PT>15: continue
        #if not abs(jet.Eta)<5: continue
        if debugmode:
            if abs(jet.PT-1069.61)<1: 
                gmethods = dir(jet)
                for ipart, part in enumerate(jet.Particles):
                    print ijet, 'genjetparticle', ipart, part.PT, part.Eta, 'pid', part.PID
                for gmethod in gmethods[:40]:
                    try:  print 'method', gmethod+':', getattr(jet, gmethod)
                    except: pass
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))
        genmet-=tlvjet
        ujet = UsefulJet(tlvjet)
        allgenjets.push_back(ujet)
        if calcMinDr(acme_objects, tlvjet, 0.3)<0.3: 
            if debugmode:
                print 'throwing away gen jet with pT, eta', jet.PT, jet.Eta, jet.Phi
                print 'gen dr was', calcMinDr(acme_objects, tlvjet, 0.01)#can check for ambiguity
                print 'while those acme objects were:'
                for obj in acme_objects: print obj.Pt(), obj.Eta(), obj.Phi()
            continue		
        genjets_.push_back(tlvjet)

    gSt = getHt(genjets_,AnHardMetJetPtCut)
    for obj in acme_objects: gSt+=obj.Pt()
    fillth1(hSt, gSt,1)

    if debugmode:
        for genthingy in genthingies:
            deltarmin = 999
            matched = False
            for genjet in allgenjets:
                deltarmin = min(deltarmin,genthingy[0].DeltaR(genjet.tlv))
                if deltarmin<0.4:
                    matched = True
                    break
            if not matched:
                if not genthingy[0].Pt()>2: continue
                print 'not matched:', genthingy[0].Pt(), genthingy[0].Eta(), genthingy[0].Phi(), 'pid =', genthingy[2], 'mindr('+str(deltarmin)+')'


    matchedCsvVec = createMatchedCsvVector(genjets_, recojets)
    genjets = CreateUsefulJetVector(genjets_, matchedCsvVec)

    ##a few global objects
    MetVec = TLorentzVector()
    MetVec.SetPtEtaPhiE(branchMissingET[0].MET,0,branchMissingET[0].Phi,branchMissingET[0].MET)

    ##observed histogram
    if debugmode:        
        print 'genjets:'
        for ijet, jet in enumerate(genjets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()
        print 'recojets:'
        for ijet, jet in enumerate(branchJets):
            print ijet, jet.PT, jet.Eta, jet.Phi
        
    tHt = getHt(recojets,AnHardMetJetPtCut,2.4)            
    tHt5 = getHt(recojets,AnHardMetJetPtCut,5.0)

    tSt = tHt + getHt(acme_objects,AnHardMetJetPtCut,2.4)    
    tSt5 = tHt5 + getHt(acme_objects,AnHardMetJetPtCut,5.0)

    if tSt>0: stratio = tSt5/tSt
    else: stratio = -1


    tHardMetVec = getHardMet(recojets,AnHardMetJetPtCut, mhtjetetacut)
    tHardMhtVec = tHardMetVec.Clone()
    tHardMetVec-=AcmeVector
    tHardMetPt, tHardMetPhi = tHardMetVec.Pt(), tHardMetVec.Phi()
    tHardMhtPt, tHardMhtPhi = tHardMhtVec.Pt(), tHardMhtVec.Phi()

    mindphi = 4
    for jet in recojets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(tHardMetVec)))

    #met_consistency = abs(branchMissingET[0].MET-tHardMetPt)/tHardMetPt
    #if not met_consistency<0.5: continue
    if abs(branchMissingET[0].MET-tHardMetPt)<100:
        if searchmode or debugmode: 
            continue        
    else:
        print ientry, 'these things are different', branchMissingET[0].MET, tHardMetPt
        #if not searchmode: continue
    

    
    #print ientry, 'these things, branchMissingET[0].MET-tHardMetPt, are different', branchMissingET[0].MET, tHardMetPt
    #if not tSt>tHardMetPt: continue



    #print 'now processing c', ientry

    if searchmode:
        print ientry, 'bogey incoming'
        print ientry, 'the many mets are, on this joyous of occasions with photons', 'acme len=',len(acme_objects)
        print 'branchMissingET[0]', branchMissingET[0].MET, 'genmet', branchGenMissingET[0].MET, 'tHardMetPt',tHardMetPt
        print 'the reco and gen st2p4s are', tSt, gSt
        print 'raw gen jet mht pt = ', getHardMet(allgenjets, 30).Pt(), 'raw hand met =', rawgenmhtthingy.Pt()
        for genjet in allgenjets:
            if not genjet.Pt()>100: continue
            drmin = 999
            matched = False
            for recojet in recojets:

                dr = recojet.DeltaR(genjet)
                if dr<drmin:
                    drmin = dr
                if drmin<0.4:
                    matched = True

            if not matched:
                for recojet in recojets:

                    dr = recojet.DeltaR(genjet)
                    print 'considering jet', dr, recojet.Pt(), recojet.Eta(), recojet.Phi()				
                    if dr<drmin:
                        drmin = dr
                    if drmin<0.4:
                        matched = True				
                print ientry, 'not matched', genjet.Pt(), genjet.Eta(), genjet.Phi()
                print 'gen mht thingy', genmet.Pt()
                exit(0)

        continue



    if tSt>0: tMetSignificance = tHardMetPt/TMath.Sqrt(tSt)
    else: tMetSignificance = 8
    tNJets = countJets(recojets,AnHardMetJetPtCut)
    tBTags = countBJets(recojets,AnHardMetJetPtCut)


    #print 'evt. number', ientry, 'nel =', len(recoelectrons), 'nmu =', len(recomuons), 'HT =', tHt


    fv = [tSt,tHt,tHardMetPt,tNJets,tBTags,int(recophotons.size()),phopt, phoeta, mindphi, stratio]


    if searchmode: print ientry, fv
    if debugmode: 
        print 'jets:'
        for ijet, jet in enumerate(recojets): print ijet, 'reco', jet.Pt(), jet.Eta(), 'phi=', jet.Phi()

        for ijet, jet in enumerate(genjets): print ijet, 'gen', jet.Pt(), jet.Eta(),  'phi=', jet.Phi()
    #if debugmode: continue	


    for regionkey in regionCuts:
        for ivar, varname in enumerate(varlist):
            hname = regionkey+'_'+varname
            if selectionFeatureVector(fv,regionkey,varname,''): 
                fillth1(histoStructDict[hname].Observed, fv[ivar], weight*genweight)

    if mktree:
        if fv[1]>=150 and fv[1]<=200 and fv[0]>500 and fv[2]>3:
            growTree(littletree, fv, jetPhis, weight)            

    #if tHt>0 and (tHardMetPt>50 and tMetSignificance>3):


fnew.cd()
writeHistoStruct(histoStructDict, 'observed')
hGenMetGenHardMetRatio.Write()
hSt.Write()
hStWeighted.Write()

hPassFit.Write()
hTotFit.Write()

hPtGenLeptons.Write()
hEtaGenLeptons.Write()
hPtRecoLeptons.Write()
hEtaRecoLeptons.Write()
hPtRecoIsoLeptons.Write()
hEtaRecoIsoLeptons.Write()
hNLepton.Write()
if mktree:
    treefile.cd()
    littletree.Write()
    treefile.Close()

print 'just created', fnew.GetName()
fnew.Close()



