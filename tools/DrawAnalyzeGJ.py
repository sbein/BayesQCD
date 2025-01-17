from ROOT import *
from utils import *
gROOT.SetBatch(1)

chain = TChain('TreeMaker2/PreSelection')
chain.Add('/eos/uscms//store/user/sbein/RebalanceAndSmear/Run2ProductionV17/*Summer16v3.GJets_DR-0p4_HT*.root')#*HT-400
chain.Show(0)
print 'nevents =', chain.GetEntries()

universalconstraint = ' abs(MET-HardMETPt)<90'
#universalconstraint = ' HardMETPt/MET<1.3'
promptname = 'Photons_nonPrompt'
#promptname = 'Photons_genMatched'
WP = 'Medium/2'
WP = 'Loose'

plotBundle = {}

plotBundle['Baseline1PhoR_HardMet'] = ['HardMETPt>>hadc(10,130,530)','NJets>1 && NPhotons'+WP+'==1 && '+promptname+'[0]==0',True]
plotBundle['Baseline1PhoR_NJets'] = ['NJets>>hadc(10,0,10)','HardMETPt>120 && NPhotons'+WP+'==1 && '+promptname+'[0]==0',True]
plotBundle['Baseline1PhoR_PhoPt'] = ['Photons[0].Pt()>>hadc(20,50,550)','HardMETPt>120 && NPhotons'+WP+'==1  && '+promptname+'[0]==0',True]

plotBundle['Baseline1PhoF_HardMet'] = ['HardMETPt>>hadc(10,130,530)','NJets>1 && NPhotons'+WP+'==1 && '+promptname+'[0]==1',True]
plotBundle['Baseline1PhoF_NJets'] = ['NJets>>hadc(10,0,10)','HardMETPt>120 && NPhotons'+WP+'==1  && '+promptname+'[0]==1',True]
plotBundle['Baseline1PhoF_PhoPt'] = ['Photons[0].Pt()>>hadc(20,50,550)','HardMETPt>120 && NPhotons'+WP+'==1  && '+promptname+'[0]==1',True]

plotBundle['Baseline1PhoR_DPhiPhoMet'] = ['min(abs(Photons[0].Phi()-HardMETPhi),2*3.14159-abs(Photons[0].Phi()-HardMETPhi))>>hadc(8,0,3.2)','HardMETPt>120 && NPhotons'+WP+'==1   && '+promptname+'[0]==0',True]
plotBundle['Baseline1PhoF_DPhiPhoMet'] = ['min(abs(Photons[0].Phi()-HardMETPhi),2*3.14159-abs(Photons[0].Phi()-HardMETPhi))>>hadc(8,0,3.2)','HardMETPt>120 && NPhotons'+WP+'==1   && '+promptname+'[0]==1',True]

'''
plotBundle['Baseline2Pho_HardMet'] = ['HardMETPt>>hadc(10,130,530)','NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_NJets'] = ['NJets>>hadc(10,0,10)','HardMETPt>120 && NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_Pho1Pt'] = ['Photons[0].Pt()>>hadc(10,50,550)','HardMETPt>120 && NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_Pho2Pt'] = ['Photons[1].Pt()>>hadc(10,50,550)','HardMETPt>120 && NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_DPhiGG'] = ['min(abs(Photons[0].Phi()-Photons[1].Phi()),2*3.14159-abs(Photons[0].Phi()-Photons[1].Phi()))>>hadc(8,0,3.2)','HardMETPt>120 && NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_DPhiPho1Met'] = ['min(abs(Photons[0].Phi()-HardMETPhi),2*3.14159-abs(Photons[0].Phi()-HardMETPhi))>>hadc(8,0,3.2)','HardMETPt>120 && NPhotons'+WP+'==2',True]
plotBundle['Baseline2Pho_DPhiPho2Met'] = ['min(abs(Photons[1].Phi()-HardMETPhi),2*3.14159-abs(Photons[1].Phi()-HardMETPhi))>>hadc(8,0,3.2)','HardMETPt>120 && NPhotons'+WP+'==2',True]
'''


fnew = TFile('ezplots_gjets.root', 'recreate')
c1 = mkcanvas()
c2 = mkcanvas('c2')

for key in plotBundle:
    drawarg, constraint, usernsvalue = plotBundle[key]
    obsweight = '1*('+constraint + ' && '+ universalconstraint + ' && IsUniqueSeed==1)'#puWeight * CrossSection
    print 'drawing', drawarg, ', with constraint:', obsweight
    chain.Draw(drawarg,obsweight)
    hobs = chain.GetHistogram().Clone(key+'_obs')
    hobs.GetYaxis().SetRangeUser(0.01,10000*hobs.GetMaximum())
    if not 'Vs' in key:
        c1.cd()
        if usernsvalue: 
            drawarg = drawarg.replace('METPt','METPtRandS').replace('NJets','NJetsRandS').replace('BTags','BTagsRandS')
            randsconstraint = constraint.replace('METPt','METPtRandS').replace('NJets','NJetsRandS').replace('BTags','BTagsRandS')
        methweight = '1/NSmearsPerEvent*('+ randsconstraint + ' && '+universalconstraint+')'#puWeight * CrossSection
        print 'drawing', drawarg, ', with constraint:', methweight
        chain.Draw(drawarg, methweight)
        hrands = chain.GetHistogram().Clone(key+'_rands') 
        hrands.GetYaxis().SetRangeUser(0.01,10000*hrands.GetMaximum())
        histoStyler(hrands, kAzure-8)
        hrands.SetFillColor(hrands.GetLineColor())
        hrands.SetFillStyle(1001)

        leg = mklegend(x1=.45, y1=.57, x2=.95, y2=.74, color=kWhite)
        hobs.SetTitle('Summer16 QCD observed')
        hobs.GetXaxis().SetTitle(key.split('_')[-1])
        hrands.GetXaxis().SetTitle(key.split('_')[-1])
        hrands.SetTitle('QCD sim. rebalance and smeared')
        hratio, hmethodsyst = FabDrawSystyRatio(c1,leg,hobs,[hrands],datamc='MC',lumi='n/a', title = '', LinearScale=False, fractionthing='truth / method')
        c1.Update()
        c1.Write('c_'+key)        
    else:
        c2.cd()
        hobs.Draw('colz')
        c2.Update()
        hobs.Write()
        c2.Write('c_'+key)        



c2 = mkcanvas('c2')

drawarg = "Jets[JetsRebalanced_origIdx[0]].Pt()/GenJets[0].Pt()>>hadc(100,0,3)"
constraint = "IsUniqueSeed==1 && fabs(JetsRebalanced[0].Eta()-GenJets[0].Eta())<0.4 && HardMETPt>100"
print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histJetResponse = chain.GetHistogram().Clone('histJetResponse')
histoStyler(histJetResponse, kGray+2)
print 'RMS', histJetResponse.GetRMS()

drawarg = "JetsRebalanced[0].Pt()/GenJets[0].Pt()>>hadc(100,0,3"
constraint = "IsUniqueSeed==1 && fabs(JetsRebalanced[0].Eta()-GenJets[0].Eta())<0.4 && HardMETPt>100"
print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histRebalanced = chain.GetHistogram().Clone('hRebalanced')
histoStyler(histRebalanced, kAzure+1)
print 'RMS', histRebalanced.GetRMS()

histJetResponse.GetXaxis().SetTitle('pT(RECO)/pT(GEN)')
histJetResponse.Draw('hist')
histRebalanced.Draw('hist same e')

leg = mklegend()
leg.AddEntry(histRebalanced, 'rebalanced jet')
leg.AddEntry(histJetResponse, 'original jet')
leg.Draw()
c2.Update()
c2.Write('response')
histJetResponse.Write()
histRebalanced.Write()

print 'just created', fnew.GetName()
fnew.Close()
