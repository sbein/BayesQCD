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
lhdHardMetJetPtCut = 15.0
AnHardMetJetPtCut = 30.0
cutoff = 15.0
isdata = False
rebalancedMetCut = 90
rebalancedMetCut = 120
#rebalancedMetCut = 150
#rebalancedMetCut = 100

#rebalancedMetCut = 180


debugmode = False
searchmode = False
usefullprior = False
usefulllikelihood = False

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
parser.add_argument("-printevery", "--printevery", type=int, default=1000,help="short run")
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
regionCuts['Baseline0Pho']      = [[0,Inf],[120,Inf],[120,Inf],[1,Inf],[0,0],  [0,0],    [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]]
regionCuts['HighHtBaseline0Pho']=[[0,Inf],[500,Inf],[120,Inf],[1,Inf],[0,0],  [0,0],    [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]]
regionCuts['Baseline0PhoLDP']      = [[0,Inf],[120,Inf],[120,Inf],[1,Inf],[0,0],  [0,0],    [-Inf,Inf],[-Inf,Inf],[0,0.5],[-Inf,Inf]]
#regionCuts['Baseline1Pho']      = [[0,Inf],[120,Inf],[120,Inf],[1,Inf],[0,0],  [1,1],    [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]]
#regionCuts['Baseline2Pho']      = [[0,Inf],[120,Inf],[120,Inf],[1,Inf],[0,0],  [2,Inf],  [-Inf,Inf],[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]]


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
templateFileName = 'usefulthings/llhd-prior-coarse-1p2width1p2tail.root'
templateFileName = 'usefulthings/llhd-prior-coarse-1p2width.root'
templateFileName = 'usefulthings/llhd-prior-coarse-2width.root'
templateFileName = 'usefulthings/llhd-prior-coarse-widenedecalhcal.root'
templateFileName = 'usefulthings/llhd-prior-coarse-widenedecalhcal_step6.root'
templateFileName = 'usefulthings/llhd-prior-coarse-widenedecalhcal_step7.root'


if usefulllikelihood:
    templateFileName = 'usefulthings/ResponseFunctionsMC17AllFilters_deepCsv.root'
    ftemplate = TFile(templateFileName)
    
else:
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
newname = 'posterior-'+infileID+'.root'
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

#GleanTemplatesFromFile(ftemplate)
print 'listing contents of prior:'
fprior.ls('splines/*')

#GleanTemplatesFromFile(ftemplate, fprior)
GleanTemplatesFromFile(fprior)




histoStructDict = {}
for region in regionCuts:
    for var in varlist:
        histname = region+'_'+var
        histoStructDict[histname] = mkHistoStruct(histname, binning)

xsec_times_lumi_over_nevents = 1.0
t0 = time.time()

for ientry in range(n2process):

    if debugmode:
        #if not ientry>122000: continue
        if not ientry in [5377, 9548, 33327, 44767]: continue
        a = 2

    c.GetEntry(ientry)
    
    if ientry%printevery==0:
        print "processing event", ientry, '/', n2process, 'time', time.time()-t0, branchMissingET[0].MET



    weight = xsec_times_lumi_over_nevents

    #if not branchMissingET[0].MET>175: continue

    acme_objects = vector('TLorentzVector')()
    recophotons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using


    ngpho=0
    if (searchmode and debugmode):
        for igp_, gp_ in enumerate(branchParticles):            
            if not gp_.E>1: continue
            if abs(gp_.PID) in [22,12,14,16,15] and gp_.PT>20:

                if abs(gp_.PID)==22: 
                    print ientry, igp_, 'got gen photon', round(gp_.PT,1), 'eta='+str(round(gp_.Eta,1)), round(gp_.Phi,2), 'pid='+str(gp_.PID), 'status='+str(gp_.Status), 'mass',round(gp_.Mass,2), 'm1', gp_.M1, 'd1', gp_.D1, 'd2', gp_.D2
                    if gp_.Status==1: ngpho+=1
                else:
                    print ientry, igp_, 'got neutrino', round(gp_.PT,1), 'eta='+str(round(gp_.Eta,1)), round(gp_.Phi,2), 'pid='+str(gp_.PID), 'status='+str(gp_.Status), 'mass',round(gp_.Mass,2), 'mother id', branchParticles[gp_.M1].PID, 'mother mother id',  branchParticles[branchParticles[gp_.M1].M1].PID

    #if not len(branchPhoton)>0: continue
    nspikes = 0
    for ipho, pho in enumerate(branchPhoton):
        if debugmode: print 'we got fresh photon', pho.PT, pho.Eta, pho.Phi
        if not pho.PT>30: continue ####
        if not abs(pho.Eta)<5.0: continue
        tlvpho = TLorentzVector()
        tlvpho.SetPtEtaPhiE(pho.PT, pho.Eta, pho.Phi, pho.E)        
        if debugmode:
            print ientry, 'acme photon', pho.PT, pho.Eta, pho.Phi, pho.IsolationVar
            #for method in dir(pho):
            #    print method+':', getattr(pho,method)
            for igp_, gp_ in enumerate(branchParticles):
                if not gp_.PT>20: continue
                if not gp_.Status==1: continue
                gp = TLorentzVector()
                gp.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
                dr = gp.DeltaR(tlvpho)
                if dr<0.1:
                    print igp_, '('+str(gp_.PID)+', status='+str(gp_.Status)+')', 'particle', gp.Pt(), gp.Eta(), gp.Phi(), 'DR=', dr

        isspike = True
        for igp_, gp_ in enumerate(branchParticles):
            if not gp_.PT>10: continue
            if not gp_.Status==1: continue
            if debugmode: print 'we got some gen particle', gp_.PT, gp_.Eta, gp_.Phi
            gp = TLorentzVector()
            gp.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
            dr = gp.DeltaR(tlvpho)
            if dr<0.1:
                isspike = False
                break
                #print igp_, '('+str(gp_.PID)+', status='+str(gp_.Status)+')', 'particle', gp.Pt(), gp.Eta(), gp.Phi(), 'DR=', dr
        if isspike: 
            nspikes+=1
            print ientry, 'tossing a spike', pho.PT, pho.Eta, 'iso:',pho.IsolationVar
            continue

        if abs(pho.Eta)<2.5 and pho.PT>50:
            recophotons.push_back(tlvpho)
            acme_objects.push_back(tlvpho)


    if not nspikes==0: continue #####spikes!
    ###if not len(recophotons)==0: continue

    if len(recophotons)>0: 
        phopt, phoeta = recophotons[0].Pt(), recophotons[0].Eta()
    else:
        phopt, phoeta = -1, -11
    '''
    if not len(branchPhoton)>0: continue
    for igp_, gp_ in enumerate(branchParticles):
        if not gp_.PT>30: continue
        if not gp_.Status==1: continue
        if not gp_.PID==22: continue
        tlvpho = TLorentzVector()
        tlvpho.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
        if debugmode:
            print ientry, 'acme photon', tlvpho.Pt(), tlvpho.Eta(), tlvpho.Phi()			
        acme_objects.push_back(tlvpho)
        if not abs(gp_.Eta)<2.4: continue		
        recophotons.push_back(tlvpho)
    '''

    if len(recophotons)>1: dphi = abs(recophotons[0].DeltaPhi(recophotons[1]))
    else: dphi = -1


    recoelectrons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for iel, el in enumerate(branchElectron):
        if not el.PT>10: continue
        #if not el.IsolationVar==0: continue
        tlvel = TLorentzVector()
        tlvel.SetPtEtaPhiE(el.PT, el.Eta, el.Phi, el.PT*TMath.CosH(el.Eta))
        if debugmode:
            print ientry, 'acme electron pt and iso', el.PT, el.IsolationVar
        if abs(el.Eta)<2.4: 
            recoelectrons.push_back(tlvel)
            acme_objects.push_back(tlvel)

    if not len(recoelectrons)==0: continue

    recomuons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for imu, mu in enumerate(branchMuon):
        if not mu.PT>10: continue
        if not abs(mu.Eta)<5.0: continue
        #if not mu.IsolationVar==0: continue
        tlvmu = TLorentzVector()
        tlvmu.SetPtEtaPhiE(mu.PT, mu.Eta, mu.Phi, mu.PT*TMath.CosH(mu.Eta))
        if debugmode:
            print ientry, 'acme muon pt and iso', mu.PT, mu.IsolationVar			
        if abs(mu.Eta)<2.5: 
            recomuons.push_back(tlvmu)
            acme_objects.push_back(tlvmu)


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
        if not jet.PT>15: continue
        if not abs(jet.Eta)<5: continue
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))
        ujet = UsefulJet(tlvjet, jet.BTag)
        if calcMinDr(acme_objects, ujet, 0.3)<0.3:
            if debugmode: print 'throwing away reco jet', ujet.Pt()
            continue
        if calcMinDr(acme_objects, ujet, 0.3)<0.55:
            onlygoodjets = False
        recojets.push_back(ujet)
        if ujet.Pt()>50 and jet.EhadOverEem>100: 
            onlygoodjets = False
            break
    if not onlygoodjets: continue



    tHardMetVec = getHardMet(recojets,AnHardMetJetPtCut, mhtjetetacut) # for debugging
    tHardMhtVec = tHardMetVec.Clone() # for debugging
    tHardMetVec-=AcmeVector # for debugging
    tHardMetPt, tHardMetPhi = tHardMetVec.Pt(), tHardMetVec.Phi() # for debugging


    if searchmode:
    	if not tHardMetPt>200: continue # for debugging        




    allgenjets = vector('UsefulJet')()
    genjets_ = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #user should interface with the input format of the event using
    for ijet, jet in enumerate(branchGenJet):
        #if not jet.PT>15: continue
        #if not abs(jet.Eta)<5: continue
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))
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
        for ijet, jet in enumerate(recojets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()            

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
    for jet in recojets:
        if not jet.Pt()>30: continue
        if not abs(jet.Eta())<2.4: continue
        mindphi = min(mindphi, abs(jet.DeltaPhi(tHardMetVec)))

    #met_consistency = abs(branchMissingET[0].MET-tHardMetPt)/tHardMetPt
    #if not met_consistency<0.5: continue

    if not abs(branchMissingET[0].MET-tHardMetPt)<100: continue
    if not tSt>tHardMetPt: continue


    ngentaus = 0
    for igp_, gp_ in enumerate(branchParticles):            
            if not gp_.PT>10: continue
            if abs(gp_.PID) ==15: 
                ngentaus+=1
                break
    if ngentaus>0: 
        print 'skipping tau'
        continue


    if tHardMetPt>120:
        print ientry, 'bogey incoming'
        print ientry, 'the many mets are, on this joyous of occasions with photons', 'acme len=',len(acme_objects)
        print 'branchMissingET[0]', branchMissingET[0].MET, 
        print 'tHardMetPt',tHardMetPt




    if tSt>0: tMetSignificance = tHardMetPt/TMath.Sqrt(tSt)
    else: tMetSignificance = 8
    tNJets = countJets(recojets,AnHardMetJetPtCut)
    tBTags = countBJets(recojets,AnHardMetJetPtCut)


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
                fillth1(histoStructDict[hname].Observed, fv[ivar], weight)

    if mktree:
        if fv[1]>=150 and fv[1]<=200 and fv[0]>500 and fv[2]>3:
            growTree(littletree, fv, jetPhis, weight)            

    fitsucceed = RebalanceJets(recojets)
    rebalancedJets = _Templates_.dynamicJets

    mHt = getHt(rebalancedJets,AnHardMetJetPtCut,2.4)
    mHt5 = getHt(rebalancedJets,AnHardMetJetPtCut,5.0)
    mSt = mHt + getHt(acme_objects,AnHardMetJetPtCut,2.4)    
    mSt5 = mHt5 + getHt(acme_objects,AnHardMetJetPtCut,5.0)


    if mSt>0: stratio = mSt5/mSt
    else: stratio = -1

    mHardMetVec = getHardMet(rebalancedJets,AnHardMetJetPtCut, mhtjetetacut)
    mHardMetVec-=AcmeVector
    mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()
    if mSt>0: mMetSignificance = mHardMetPt/TMath.Sqrt(mSt)
    else: mMetSignificance = 8	

    mNJets = countJets(rebalancedJets,AnHardMetJetPtCut)
    mBTags = countBJets(rebalancedJets,AnHardMetJetPtCut)###

    hope = (fitsucceed and mHardMetPt<rebalancedMetCut)# mHardMetPt>min(mSt/2,180):# was 160	

    redoneMET = redoMET(MetVec,recojets,rebalancedJets)
    mMetPt,mMetPhi = redoneMET.Pt(), redoneMET.Phi()

    mindphi = 4
    for jet in rebalancedJets: 
        if not jet.Pt()>30: continue
        if not abs(jet.Eta())<2.4: continue
        mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))

    fv = [mSt,mHt,mHardMetPt,mNJets,mBTags, int(recophotons.size()),phopt, phoeta, mindphi, stratio]

    for regionkey in regionCuts:      
        for ivar, varname in enumerate(varlist):
            hname = regionkey+'_'+varname
            if selectionFeatureVector(fv,regionkey,varname,''): 
                fillth1(histoStructDict[hname].Rebalanced, fv[ivar],weight)

    fillth1(hTotFit, fv[3], weight)

    if hope: fillth1(hPassFit, fv[3], weight)

    nsmears = 20*bootupfactor
    weight = xsec_times_lumi_over_nevents / nsmears


    for i in range(nsmears):

        if not hope: 
            break
        RplusSJets = smearJets(rebalancedJets,99+_Templates_.nparams)
        rpsHt = getHt(RplusSJets,AnHardMetJetPtCut,2.4)
        rpsHt5 = getHt(RplusSJets,AnHardMetJetPtCut,5.0)

        rpsSt = rpsHt + getHt(acme_objects,AnHardMetJetPtCut,2.4)
        rpsSt5 = rpsHt5 + getHt(acme_objects,AnHardMetJetPtCut,5.0)

        if rpsSt>0: stratio = rpsSt5/rpsSt
        else: stratio = -1


        rpsHardMetVec = getHardMet(RplusSJets,AnHardMetJetPtCut, mhtjetetacut)
        rpsHardMetVec-=AcmeVector
        rpsHardMetPt, rpsHardMetPhi = rpsHardMetVec.Pt(), rpsHardMetVec.Phi()
        if rpsSt>0: rpsMetSignificance = rpsHardMetPt/TMath.Sqrt(rpsSt)
        else: rpsMetSignificance = 8			
        rpsNJets = countJets(RplusSJets,AnHardMetJetPtCut)
        rpsBTags = countBJets(RplusSJets,AnHardMetJetPtCut)
        mindphi = 4
        for jet in RplusSJets: 
            if not jet.Pt()>30: continue
            if not abs(jet.Eta())<2.4: continue            
            mindphi = min(mindphi, abs(jet.DeltaPhi(rpsHardMetVec)))            
        fv = [rpsSt,rpsHt,rpsHardMetPt,rpsNJets,rpsBTags, int(recophotons.size()),phopt, phoeta, mindphi, stratio]
        for regionkey in regionCuts:     
            for ivar, varname in enumerate(varlist):
                hname = regionkey+'_'+varname
                if selectionFeatureVector(fv,regionkey,varname,''):
                    fillth1(histoStructDict[hname].RplusS, fv[ivar],weight)
        if mktree and 'JetHT' in physicsProcess and fv[1]>=150 and fv[1]<200 and fv[0]>=500 and fv[2]>=4: 
            growTree(littletree, fv, jetPhis, weight)


    if isdata: continue    
    genMetVec = mkmet(branchGenMissingET[0].MET, branchGenMissingET[0].Phi)
    weight = xsec_times_lumi_over_nevents

    gHt = getHt(genjets,AnHardMetJetPtCut,2.4)
    gHt5 = getHt(genjets,AnHardMetJetPtCut,5.0)    
    gSt = gHt + getHt(acme_objects,AnHardMetJetPtCut,2.4)
    gSt5 = gHt5 + getHt(acme_objects,AnHardMetJetPtCut,5.0)

    if gSt>0: stratio = gSt5/gSt
    else: stratio = -1

    gHardMetVec = getHardMet(genjets,AnHardMetJetPtCut, mhtjetetacut)
    gHardMetVec-=AcmeVector
    gHardMetPt, gHardMetPhi = gHardMetVec.Pt(), gHardMetVec.Phi()
    if gSt>0: gMetSignificance = gHardMetPt/TMath.Sqrt(gSt)
    else: gMetSignificance = 8	


    if searchmode:
        print 'gHardMetPt', gHardMetPt
        print 'branchGenMissingET', branchGenMissingET[0].MET

        gMetFromGenJets = TLorentzVector()
        for gj in allgenjets: gMetFromGenJets-=gj.tlv

        gUnadulturatedGenHardMetVec = getHardMet(allgenjets,AnHardMetJetPtCut, mhtjetetacut)

        print 'gMetFromAllGenJets', gMetFromGenJets.Pt()
        print 'gUnadulturatedGenHardMetVec.Pt()', gUnadulturatedGenHardMetVec.Pt()###ready to try this

        if debugmode:
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


    ###Delphes filter, because I think Delphes is mis-computing its own gen MET
    if gHardMetPt>80 and False: 
        fillth1(hGenMetGenHardMetRatio,abs(gHardMetPt-genMetVec.Pt())/gHardMetPt)
        if abs(gHardMetPt-genMetVec.Pt())/gHardMetPt>0.5: 
            print 'skipping', ientry, gHardMetPt, genMetVec.Pt()
            #continue
    ###Delphes

    gNJets = countJets(genjets,AnHardMetJetPtCut)
    gBTags = countBJets(genjets,AnHardMetJetPtCut)


    if debugmode and len(acme_objects)>0:
        print 'going in there'
        #c.Show(ientry)
        print 'ientry', ientry, [tSt,tHardMetPt,tNJets,tBTags, int(recophotons.size()),phopt, phoeta]
        print 'acme', acme_objects.size()
        print 'missingET', branchMissingET[0].MET
        print 'gen missingET', genMetVec.Pt(), genMetVec.Phi()
        print 'genHardMet', gHardMetVec.Pt(), gHardMetVec.Phi()	
        hmmod1 = gHardMetVec.Clone()
        hmmod1-=acme_objects[0]
        print 'hmmod1', hmmod1.Pt(), hmmod1.Phi()		
        hmmod2 = gHardMetVec.Clone()
        hmmod2+=acme_objects[0]
        print 'hmmod2', hmmod2.Pt(), hmmod2.Phi()				
        print ''
        print ''	
        print 'genjets'
        for ijet, jet in enumerate(genjets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()

        print 'recoHardMet', tHardMetVec.Pt(), tHardMetVec.Phi()
        print 'endgame recojets:'
        for ijet, jet in enumerate(recojets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()
            matchedgen = calcMinDr(genjets, jet, 0.01)
            print 'matchedgenjet pt', matchedgen.Pt(), matchedgen.Eta(), matchedgen.Phi(), 'dr', matchedgen.DeltaR(jet.tlv)
            print 'particles in der naher:'
            for igp_, gp_ in enumerate(branchParticles):
                break
                #if not gp_.PT>10: continue
                if not gp_.Status==1: continue
                gp = TLorentzVector()
                gp.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
                dr = gp.DeltaR(jet.tlv)
                if dr<0.8:
                    print igp_, '('+str(gp_.PID)+', status='+str(gp_.Status)+')', 'particle', gp.Pt(), gp.Eta(), gp.Phi(), 'DR=', dr

        print 'acme_objects'
        for ijet, jet in enumerate(acme_objects):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()
        '''	
        print 'other particles:'
        for ijet, jet in enumerate(branchParticles):
            part = TLorentzVector()
            if not jet.PT>10: continue
            print ijet, jet.PT, jet.Eta, jet.Phi, '('+str(jet.PID)+')', '-'+str(jet.Status)
        '''
    mindphi = 4
    for jet in genjets: mindphi = min(mindphi, abs(jet.DeltaPhi(gHardMetVec)))
    fv = [gSt,gHt,gHardMetPt,gNJets,gBTags, int(recophotons.size()),phopt, phoeta, mindphi, stratio]	    
    for regionkey in regionCuts:
        for ivar, varname in enumerate(varlist):
            hname = regionkey+'_'+varname
            if selectionFeatureVector(fv,regionkey,varname,''): 
                fillth1(histoStructDict[hname].Gen, fv[ivar],weight)

    #Gen-smearing
    nsmears = 3
    if isdata: weight = 1.0*prescaleweight/nsmears    
    else: weight = xsec_times_lumi_over_nevents / nsmears
    for i in range(nsmears):
        if not (gHardMetPt<rebalancedMetCut): break
        smearedJets = smearJets(genjets,9999)
        #smearedJets,csvSmearedJets = smearJets(genjets,matchedCsvVec,_Templates_.ResponseFunctions,_Templates_.hEtaTemplate,_Templates_.hPtTemplate,999)
        mHt = getHt(smearedJets,AnHardMetJetPtCut,2.4)
        mHt5 = getHt(smearedJets,AnHardMetJetPtCut,5.0)
        mSt = mHt + getHt(acme_objects,AnHardMetJetPtCut,2.4)
        mSt5 = mHt5 + getHt(acme_objects,AnHardMetJetPtCut,5.0)        
        if mSt>0: stratio = mSt5/mSt
        else: stratio = -1

        mHardMetVec = getHardMet(smearedJets,AnHardMetJetPtCut, mhtjetetacut)
        mHardMetVec-=AcmeVector
        mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()
        if mSt>0: mMetSignificance = mHardMetPt/TMath.Sqrt(mSt)
        else: mMetSignificance = 8		
        mNJets = countJets(smearedJets,AnHardMetJetPtCut)
        mBTags = countBJets(smearedJets,AnHardMetJetPtCut)
        redoneMET = redoMET(genMetVec, genjets, smearedJets)
        mMetPt, mMetPhi = redoneMET.Pt(), redoneMET.Phi()
        mindphi = 4
        for jet in smearedJets: 
            if not jet.Pt()>30: continue
            if not abs(jet.Eta())<2.4: continue
            mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))	
        fv = [mSt,mHt,mHardMetPt,mNJets,mBTags, int(recophotons.size()),phopt, phoeta, mindphi, stratio]

        for regionkey in regionCuts:
            for ivar, varname in enumerate(varlist):
                hname = regionkey+'_'+varname
                if selectionFeatureVector(fv,regionkey,varname,''): 
                    fillth1(histoStructDict[hname].GenSmeared, fv[ivar],weight)

fnew.cd()
writeHistoStruct(histoStructDict)
hGenMetGenHardMetRatio.Write()
hSt.Write()
hStWeighted.Write()

hPassFit.Write()
hTotFit.Write()
if mktree:
    treefile.cd()
    littletree.Write()
    treefile.Close()
print 'just created', fnew.GetName()
fnew.Close()



