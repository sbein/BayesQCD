#Welcome to the industrial age of Sam's rebalance and smear code. You're going to have a lot of fun!
import os,sys
from ROOT import *
from array import array
from glob import glob
import numpy as np
import time

#python tools/MaximizePosteriorMakeTree.py --fnamekeyword /nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/SampleProduction/delphes/rootfiles_widenedhcalStep7/delphes_Weak_98of499.root

### a few parameters to set that define the selection associated with the hard MET
hardmet_jetetacut = 5.0 # eta acceptance for jets in hard MET
hardmet_jetptcut = 30.0 # pT threshold for jets going in to hard MET
isdata = False # in case this is adapted to run over data
rebalancedMetCut = 120 #somewhat tunable maximum rebalanced hard MET value accepted as a seed (should be less than the analysis baseline)
hardMetCutForSkim = 120

##load in delphes libraries to access input collections:
gSystem.Load("/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/SampleProduction/delphes/build/libDelphes")

##load in UsefulJet class,which the rebalance and smear code uses
gROOT.ProcessLine(open('src/UsefulJet.cc').read())
exec('from ROOT import *')
gROOT.ProcessLine(open('src/BayesRandS.cc').read())
exec('from ROOT import *')

##read in command line arguments
defaultInfile_ = "/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src//SampleProduction/delphes/rootfiles_widenedhcalStep7/delphes_qcd_12of80.root"
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v","--verbosity",type=int,default=1,help="analyzer script to batch")
parser.add_argument("-printevery","--printevery",type=int,default=1000,help="short run")
parser.add_argument("-fin","--fnamekeyword",type=str,default=defaultInfile_,help="file")
parser.add_argument("-quickrun","--quickrun",type=bool,default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
verbosity = args.verbosity
printevery = args.printevery
quickrun = args.quickrun

#Delphes b-tagging binary
BTag_Cut = 0.5 

#This should be set to False if there are leptons in the analysis
vetoleptons = True

##declare and load a tree
c = TChain('Delphes')
for fname in inputFiles: 
    print ('adding '+fname)
    c.Add(fname)
nentries = c.GetEntries()
c.Show(0)
if quickrun: n2process = min(5000,nentries)
else: n2process = nentries
    
print ('n(entries) = '+str(n2process))

#feed the tree to delphes,set up which branches need to be used
treeReader = ExRootTreeReader(c)
numberOfEntries = treeReader.GetEntries()
branchHT = treeReader.UseBranch("ScalarHT")
branchJets = treeReader.UseBranch("Jet")
branchMissingET = treeReader.UseBranch("MissingET")
branchGenMissingET = treeReader.UseBranch("GenMissingET")
branchPhoton = treeReader.UseBranch("Photon")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")
branchParticles = treeReader.UseBranch("Particle")


#pick up the jet smearing templates. A default set is provided but these are experiment specific
#to templates can be remade by running 
'''
python tools/LlhdPriorHistmaker.py <input file [usually large MC QCD file or many files]>
python tools/articulateSplines.py "<path-to-output/*.root>"
'''
templateFileName = 'usefulthings/llhd-prior-coarse-realisticHcal.root'
ftemplate = TFile(templateFileName)
print ('using templates from '+ templateFileName)
hPtTemplate = ftemplate.Get('hPtTemplate')
templatePtAxis = hPtTemplate.GetXaxis()
hEtaTemplate = ftemplate.Get('hEtaTemplate')
templateEtaAxis = hEtaTemplate.GetXaxis()
fprior = TFile(templateFileName)
hHtTemplate = fprior.Get('hHtTemplate')
templateStAxis = hHtTemplate.GetXaxis()


##Create output file
infileID = fnamekeyword.split('/')[-1].replace('.root','')
newname = 'littletree_fin'+infileID.replace('*','')+'.root'
fnew = TFile(newname,'recreate')
print ('creating '+newname)

CrossSectionPb = np.zeros(1,dtype=float)

#set cross section based on the file name.
isqcd = False
isgjet = False
if '_qcd_' in fnamekeyword: 
    isqcd = True
    CrossSectionPb[0] = 7.178e-05 * 1000000000
    nsmears = 100
elif 'WZ' in fnamekeyword: 
    CrossSectionPb[0] = 2.496e-04 * 1000000000
    nsmears = 10
elif 'Top' in fnamekeyword:
    CrossSectionPb[0] = 1.811e-07 * 1000000000
    nsmears = 10
elif 'T2qq' in fnamekeyword:
    CrossSectionPb[0] = 4*0.0323744#this is for 14 TeV T2qq-1150GeV, x4 for squark degeneracy 
    nsmears = 10   
elif 'Weak' in fnamekeyword:
    CrossSectionPb[0] = 1.430e-05  * 1000000000
    nsmears = 1
elif 'GJet' in fnamekeyword:
    CrossSectionPb[0] = 8.129e-07 * 1000000000
    nsmears = 100
    isgjet = True
else:
    CrossSectionPb[0] = 1
    nsmears = 0              

print ('CrossSectionPb[0] '+str(CrossSectionPb[0]))
littletree = TTree('littletree','littletree')
tcounter = TTree('tcounter','tcounter')
HT = np.zeros(1,dtype=float)
HardMetPt = np.zeros(1,dtype=float)
MinDPhiJetsHardMet = np.zeros(1,dtype=float)
Jet1Pt = np.zeros(1,dtype=float)
Jet2Pt = np.zeros(1,dtype=float)
Jet3Pt = np.zeros(1,dtype=float)
Jet4Pt = np.zeros(1,dtype=float)
Jet1Eta = np.zeros(1,dtype=float)
Jet2Eta = np.zeros(1,dtype=float)
Jet3Eta = np.zeros(1,dtype=float)
Jet4Eta = np.zeros(1,dtype=float)
Jet1DPhi = np.zeros(1,dtype=float)
Jet2DPhi = np.zeros(1,dtype=float)
Jet3DPhi = np.zeros(1,dtype=float)
Jet4DPhi = np.zeros(1,dtype=float)
NJets = np.zeros(1,dtype=int)
BTags = np.zeros(1,dtype=int)
NSmearsPerEvent = np.zeros(1,dtype=int)
IsRandS = np.zeros(1,dtype=int)
NLeps = np.zeros(1,dtype=int)    
FitSucceed = np.zeros(1,dtype=int)    

littletree.Branch('CrossSectionPb',CrossSectionPb,'CrossSectionPb/D')    
Jets = std.vector('TLorentzVector')()
littletree.Branch('Jets',Jets)
jetsRebalanced = std.vector('TLorentzVector')()
littletree.Branch('JetsRebalanced',jetsRebalanced)
Jets_btag = std.vector('int')()
littletree.Branch('Jets_btag',Jets_btag)
littletree.Branch('HardMetPt',HardMetPt,'HardMetPt/D')
littletree.Branch('HT',HT,'HT/D')
littletree.Branch('MinDPhiJetsHardMet',MinDPhiJetsHardMet,'MinDPhiJetsHardMet/D')
littletree.Branch('Jet1Pt',Jet1Pt,'Jet1Pt/D')
littletree.Branch('Jet2Pt',Jet2Pt,'Jet2Pt/D')
littletree.Branch('Jet3Pt',Jet3Pt,'Jet3Pt/D')
littletree.Branch('Jet4Pt',Jet4Pt,'Jet4Pt/D')
littletree.Branch('Jet1DPhi',Jet1DPhi,'Jet1DPhi/D')
littletree.Branch('Jet2DPhi',Jet2DPhi,'Jet2DPhi/D')
littletree.Branch('Jet3DPhi',Jet3DPhi,'Jet3DPhi/D')
littletree.Branch('Jet4DPhi',Jet4DPhi,'Jet4DPhi/D')
littletree.Branch('Jet1Eta',Jet1Eta,'Jet1Eta/D')
littletree.Branch('Jet2Eta',Jet2Eta,'Jet2Eta/D')
littletree.Branch('Jet3Eta',Jet3Eta,'Jet3Eta/D')
littletree.Branch('Jet4Eta',Jet4Eta,'Jet4Eta/D')
littletree.Branch('BTags',BTags,'BTags/I')
littletree.Branch('NJets',NJets,'NJets/I')
littletree.Branch('NSmearsPerEvent',NSmearsPerEvent,'NSmearsPerEvent/I')
littletree.Branch('FitSucceed',FitSucceed,'FitSucceed/I')   
littletree.Branch('IsRandS',IsRandS,'IsRandS/I') 
littletree.Branch('NLeps',NLeps,'NLeps/I')

photons = std.vector('TLorentzVector')()
littletree.Branch('photons',photons)    

#prepareLittleTree(littletree)
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
print ('listing contents of prior:')
fprior.ls('splines/*')

#GleanTemplatesFromFile(ftemplate,fprior)
GleanTemplatesFromFile(fprior)


t0 = time.time()

for ientry in range(n2process):

    c.GetEntry(ientry)

    if ientry%printevery==0:
        print ("processing event %d/%d; time remaining: %f" % (ientry,n2process,time.time()-t0))

    tcounter.Fill()
    weight = CrossSectionPb

    #the acme objects are the isolated particles in the event who's energy is presumed to be measured perfectly: ele, mu, pho, etc...
    #we don't rebalance and/or smear these objects, so they need to be kept track of independently
    acme_objects = vector('TLorentzVector')()

    #this is where the code must interface with the input format (jets collections, etc)
    
    #build up the vector of photons using TLorentzVectors;     
    recophotons = vector('TLorentzVector')()
    nspikes = 0
    for ipho,pho in enumerate(branchPhoton):
        if not pho.PT>100: continue ####
        if not abs(pho.Eta)<5.0: continue
        tlvpho = TLorentzVector()
        tlvpho.SetPtEtaPhiE(pho.PT,pho.Eta,pho.Phi,pho.E)

        #found delphes had a few unexplainable ECal spikes - these must be dealt with on an case-by-case basis
        
        isspike = False #trying to do away with this for the moment
        if isgjet:
            isspike = True
            for igp_,gp_ in enumerate(branchParticles):
                if not gp_.PT>10: continue
                if not gp_.Status==1: continue
                gp = TLorentzVector()
                gp.SetPtEtaPhiE(gp_.PT,gp_.Eta,gp_.Phi,gp_.E)
                dr = gp.DeltaR(tlvpho)
                if dr<0.1:
                    isspike = False
                    break
            if isspike: 
                nspikes+=1
                break
        if abs(pho.Eta)<2.5 and pho.PT>50:
            recophotons.push_back(tlvpho)
            acme_objects.push_back(tlvpho)

    #for now let's treat these spikes like something that can be filtered
    if not nspikes==0: continue ##spikes!

    #tau background must be measured in some other way, since it has real MET.
    if isqcd or isgjet:
        ngentaus = 0
        for igp_,gp_ in enumerate(branchParticles):            
            if not gp_.PT>10: continue
            if abs(gp_.PID) ==15: 
                ngentaus+=1
                break
        if ngentaus>0: 
            #events with taus have real MET and must be predicted by another method
            continue
        
    if len(recophotons)>0: phopt,phoeta = recophotons[0].Pt(),recophotons[0].Eta()
    else: phopt,phoeta = -1,-11

    #fill container for tree branch
    photons.clear()
    for photon in recophotons: photons.push_back(photon)

    #build vector of electrons using TLorentzVectors;         
    recoelectrons = vector('TLorentzVector')()
    for iel,el in enumerate(branchElectron):
        if not el.PT>10: continue
        #if not el.IsolationVar==0: continue
        tlvel = TLorentzVector()
        tlvel.SetPtEtaPhiE(el.PT,el.Eta,el.Phi,el.PT*TMath.CosH(el.Eta))

    #electron veto
    if vetoleptons:
        if not len(recoelectrons)==0: continue

    recomuons = vector('TLorentzVector')()
    #build vector of muons using TLorentzVectors; 
    for imu,mu in enumerate(branchMuon):
        if not mu.PT>10: continue
        if not abs(mu.Eta)<5.0: continue
        #if not mu.IsolationVar==0: continue
        tlvmu = TLorentzVector()
        tlvmu.SetPtEtaPhiE(mu.PT,mu.Eta,mu.Phi,mu.PT*TMath.CosH(mu.Eta))
        
        if abs(mu.Eta)<2.5: 
            recomuons.push_back(tlvmu)
            acme_objects.push_back(tlvmu)

    if vetoleptons:
        if not len(recomuons)==0: continue

    #to be the vectorial sum of presumed perfectly measured energy
    AcmeVector = TLorentzVector()
    AcmeVector.SetPxPyPzE(0,0,0,0)
    for obj in acme_objects: AcmeVector+=obj		

    _Templates_.AcmeVector = AcmeVector	

    ##declare empty vector of UsefulJets (in c++,std::vector<UsefulJet>):
    recojets = vector('UsefulJet')()

    #build up the vector of jets using TLorentzVectors; 
    onlygoodjets = True
    for ijet,jet in enumerate(branchJets):
        if not jet.PT>15: continue
        if not abs(jet.Eta)<5: continue
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.PT,jet.Eta,jet.Phi,jet.PT*TMath.CosH(jet.Eta))
        ujet = UsefulJet(tlvjet,jet.BTag)
        if calcMinDr(acme_objects,ujet,0.3)<0.3:
            continue
        if calcMinDr(acme_objects,ujet,0.3)<0.55:
            onlygoodjets = False
        recojets.push_back(ujet)
        if ujet.Pt()>50 and jet.EhadOverEem>100: 
            onlygoodjets = False
            break
    if not onlygoodjets: continue


    tHardMetVec = getHardMet(recojets,hardmet_jetptcut,hardmet_jetetacut) # for debugging
    tHardMhtVec = tHardMetVec.Clone() # for debugging
    tHardMetVec-=AcmeVector # for debugging
    tHardMetPt,tHardMetPhi = tHardMetVec.Pt(),tHardMetVec.Phi() # for debugging

    tmindphi = 4
    for jet in recojets:
        if not jet.Pt()>30: continue
        if not abs(jet.Eta())<2.4: continue
        tmindphi = min(tmindphi,abs(jet.DeltaPhi(tHardMetVec)))
        

    ##a few global objects
    MetVec = TLorentzVector()
    MetVec.SetPtEtaPhiE(branchMissingET[0].MET,0,branchMissingET[0].Phi,branchMissingET[0].MET)

    #H_T, scalar sum or jet pT
    tHt = getHt(recojets,hardmet_jetptcut,2.4)
    
    #S_T, scalar sum of everything
    tSt = tHt + getHt(acme_objects,hardmet_jetptcut,2.4)    


    #event filter against events with lots of low-pT activity, safe for typical DM searches
    if not abs(branchMissingET[0].MET-tHardMetPt)<100: continue
    if not tSt>tHardMetPt: continue

    #
    tNJets = countJets(recojets,hardmet_jetptcut)
    tBTags = countBJets(recojets,2.4)


    #do reblancing, which knows about global variable _Templates_
    fitsucceed = RebalanceJets(recojets)
    rebalancedJets = _Templates_.dynamicJets

    #derive any event quantities from the rebalanced jets
    mHardMetVec = getHardMet(rebalancedJets,hardmet_jetptcut,hardmet_jetetacut)
    mHardMetVec-=AcmeVector
    mHardMetPt,mHardMetPhi = mHardMetVec.Pt(),mHardMetVec.Phi()

    
    eventWasSufficientlyRebalanced = (fitsucceed and mHardMetPt<rebalancedMetCut)

    #keep track of numbers of insufficiently rebalanced events, potential weight or for systematic studies?
    hTotFit.Fill(tNJets,weight)
    if eventWasSufficientlyRebalanced: hPassFit.Fill(tNJets,weight)
        
    weight = CrossSectionPb / nsmears

    #fill the tree with unspoiled event characteristics
    if tHardMetPt>hardMetCutForSkim:
        Jets.clear()
        jetsRebalanced.clear()
        Jets_btag.clear()

        for ujet in recojets:
            Jets.push_back(ujet.tlv)

        for ujet in rebalancedJets:
            jetsRebalanced.push_back(ujet.tlv)                

        HardMetPt[0] = tHardMetPt
        MinDPhiJetsHardMet[0] = tmindphi
        HT[0] = tHt
        NJets[0] = tNJets

        if NJets[0]>0: Jet1Pt[0],Jet1Eta[0],Jet1DPhi[0] = Jets[0].Pt(),Jets[0].Eta(),Jets[0].DeltaPhi(tHardMetVec)
        else: Jet1Pt[0],Jet1Eta[0],Jet1DPhi[0] = 0,-4,-4
        if NJets[0]>1: Jet2Pt[0],Jet2Eta[0],Jet2DPhi[0] = Jets[1].Pt(),Jets[1].Eta(),Jets[1].DeltaPhi(tHardMetVec)
        else: Jet2Pt[0],Jet2Eta[0],Jet2DPhi[0] = 0,-4,-4       
        if NJets[0]>2: Jet3Pt[0],Jet3Eta[0],Jet3DPhi[0] = Jets[2].Pt(),Jets[2].Eta(),Jets[2].DeltaPhi(tHardMetVec)
        else: Jet3Pt[0],Jet3Eta[0],Jet3DPhi[0] = 0,-4,-4
        if NJets[0]>3: Jet4Pt[0],Jet4Eta[0],Jet4DPhi[0] = Jets[3].Pt(),Jets[3].Eta(),Jets[3].DeltaPhi(tHardMetVec)
        else: Jet4Pt[0],Jet4Eta[0],Jet4DPhi[0] = 0,-4,-4  

        BTags[0] = tBTags
        NSmearsPerEvent[0] = nsmears
        IsRandS[0] = 0
        NLeps[0] = len(recomuons)+len(recoelectrons)
        FitSucceed[0] = fitsucceed

        littletree.Fill()

    #smear the seed events to make "new reco data"
    for i in range(nsmears):
        if not eventWasSufficientlyRebalanced:  break
        RplusSJets = smearJets(rebalancedJets,99+_Templates_.nparams)
        rpsHt = getHt(RplusSJets,hardmet_jetptcut,2.4)
        rpsHt5 = getHt(RplusSJets,hardmet_jetptcut,5.0)

        rpsSt = rpsHt + getHt(acme_objects,hardmet_jetptcut,2.4)
        rpsSt5 = rpsHt5 + getHt(acme_objects,hardmet_jetptcut,5.0)


        rpsHardMetVec = getHardMet(RplusSJets,hardmet_jetptcut,hardmet_jetetacut)
        rpsHardMetVec-=AcmeVector
        rpsHardMetPt,rpsHardMetPhi = rpsHardMetVec.Pt(),rpsHardMetVec.Phi()
        if rpsSt>0: rpsMetSignificance = rpsHardMetPt/TMath.Sqrt(rpsSt)
        else: rpsMetSignificance = 8			
        rpsNJets = countJets(RplusSJets,hardmet_jetptcut)
        rpsBTags = countBJets(RplusSJets,2.4)
        rpsmindphi = 4
        for jet in RplusSJets: 
            if not jet.Pt()>30: continue
            if not abs(jet.Eta())<2.4: continue            
            rpsmindphi = min(rpsmindphi,abs(jet.DeltaPhi(rpsHardMetVec)))            

        #fill the tree with R&S collections
        if rpsHardMetPt>hardMetCutForSkim:
            Jets.clear()
            jetsRebalanced.clear()
            Jets_btag.clear()

            for ujet in RplusSJets:
                Jets.push_back(ujet.tlv)

            for ujet in rebalancedJets:
                jetsRebalanced.push_back(ujet.tlv)                
            HardMetPt[0] = rpsHardMetPt
            MinDPhiJetsHardMet[0] = rpsmindphi
            HT[0] = rpsHt
            NJets[0] = rpsNJets
            if NJets[0]>0: Jet1Pt[0],Jet1Eta[0],Jet1DPhi[0] = Jets[0].Pt(),Jets[0].Eta(),Jets[0].DeltaPhi(tHardMetVec)
            else: Jet1Pt[0],Jet1Eta[0],Jet1DPhi[0] = 0,-4,-4
            if NJets[0]>1: Jet2Pt[0],Jet2Eta[0],Jet2DPhi[0] = Jets[1].Pt(),Jets[1].Eta(),Jets[1].DeltaPhi(tHardMetVec)
            else: Jet2Pt[0],Jet2Eta[0],Jet2DPhi[0] = 0,-4,-4       
            if NJets[0]>2: Jet3Pt[0],Jet3Eta[0],Jet3DPhi[0] = Jets[2].Pt(),Jets[2].Eta(),Jets[2].DeltaPhi(tHardMetVec)
            else: Jet3Pt[0],Jet3Eta[0],Jet3DPhi[0] = 0,-4,-4
            if NJets[0]>3: Jet4Pt[0],Jet4Eta[0],Jet4DPhi[0] = Jets[3].Pt(),Jets[3].Eta(),Jets[3].DeltaPhi(tHardMetVec)
            else: Jet4Pt[0],Jet4Eta[0],Jet4DPhi[0] = 0,-4,-4                 
            BTags[0] = rpsBTags
            NSmearsPerEvent[0] = nsmears
            IsRandS[0] = 1
            NLeps[0] = len(recomuons)+len(recoelectrons)
            FitSucceed[0] = fitsucceed                
            littletree.Fill()
	

#write any objects needed
fnew.cd()
hGenMetGenHardMetRatio.Write()
hSt.Write()
hStWeighted.Write()

hPassFit.Write()
hTotFit.Write()

littletree.Write()
tcounter.Write()
print ('just created '+fnew.GetName())
fnew.Close()



